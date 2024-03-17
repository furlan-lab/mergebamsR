extern crate bam;
extern crate csv;
extern crate flate2;
extern crate itertools;

use itertools::Itertools;
use differ::{Differ, Tag};
use bam::RecordWriter;
use bam::record::tags::TagValue;
use std::str;
use std::error::Error;


#[derive(Clone)]
struct Params<'a> {
    bams: Vec<&'a str>,
    out_path: &'a str,
    names: Vec<Option<Vec<String>>>,
    prefixes: Vec<&'a str>,
    threads: usize,
}


pub fn mergebams_rust<'a> (bams: Vec<&str>, out_path: &'a str, names: Vec<Option<Vec<String>>>, prefixes: Vec<&str>){
    let params = load_params(bams, names, prefixes, out_path);
    let _header_result = checkheaders(params.clone());
    if let Ok((header, params)) = checkheaders(params){
        let _params = addtags(params, header);
    }else{
        eprintln!("ERROR: BAM header sequences do not match - you will need to fix this before merging bams");
    }

}


fn load_params<'a>(bams: Vec<&'a str>, names: Vec<Option<Vec<String>>>, prefixes: Vec<&'a str>, out_path: &'a str) -> Params<'a> {
    let threads = 1;
    
    Params{
        bams: bams,
        out_path: out_path,
        names: names,
        prefixes: prefixes,
        threads: threads,
    }
}


fn checkheaders(params: Params) -> Result<(bam::Header, Params), &'static str>{
    let bam_it = params.bams.iter().combinations(2);
    let mut grand_count = 0;
    for bam_pairs in bam_it {
        let mut header_vec = Vec::new();
        for inbam in bam_pairs.iter() {
            let hreader = bam::BamReader::from_path(inbam, 0).unwrap();
            let header = hreader.header().clone();
            header_vec.push(header);
        }
        let mut count = 0;
        let a = header_vec[0].reference_names();
        let b = header_vec[1].reference_names();
        let differ = Differ::new(&a, &b);
        for span in differ.spans() {
            match span.tag {
                Tag::Equal => (),
                _ => count+=1,
            }
        }

        eprintln!("Found {} discrepencies between sequence names in the header of:\n{}\nand:\n{}\n", count, bam_pairs[0], bam_pairs[1]);
        grand_count+=count;
    }
    if grand_count == 0{
        let new_header = make_new_header(params.bams.clone());
        return Ok((new_header, params));
    } else {
        return Err("Discrepent headers")
    }
    
}

fn make_new_header(bam_vec: Vec<&str>) -> bam::Header {
    // assumes no discrepencies in the headers across bams in bam_vec
    let hreader = bam::BamReader::from_path(bam_vec[0], 0).unwrap();
    let mut out_header = hreader.header().clone();
    let mergebam_line = bam_vec.join(", ");
    let _msg = out_header.push_line(&("@CO\tmergebams has included the BAM records from the following files (using the header from the first): ".to_owned()+&mergebam_line));
    return out_header;
}

fn addtags(params: Params, header: bam::Header) -> Params{
    let out_path = params.out_path.to_string()+"/out_path.bam";
    let fail_bam = params.out_path.to_string()+"/fail_bam.bam";
    let out_path_msg = out_path.clone();
    let bam_vec = params.bams.clone();
    let bam_vec_msg = params.bams.join(" and ");
    let lab_vec = params.prefixes.clone();
    let (read_threads, write_threads) = if (*&params.threads as i8) > 2{
        (((*&params.threads/2) -1) as u16, ((*&params.threads/2) -1) as u16)
    } else {
        (0 as u16, 0 as u16)
    };
    let mut fail_count = 0;
    let mut pass_count = 0;
    let mut other_count = 0;
    let mut pass_writer = bam::BamWriter::build()
        .write_header(true)
        .additional_threads(write_threads)
        .from_path(out_path, header.clone()).unwrap();
    let mut fail_writer = bam::BamWriter::build()
        .write_header(true)
        .additional_threads(0)
        .from_path(fail_bam, header.clone()).unwrap();
    eprintln!("Headers ok\nWriting:\n{}\nfrom:\n{}\n", out_path_msg, bam_vec_msg);
    for (pos, inbam) in bam_vec.iter().enumerate() {
        let reader = bam::BamReader::from_path(inbam.to_string(), read_threads).unwrap();
        let names = &params.names[pos];
        let mut filter: bool = false;
        if names.is_some(){
            filter = true;
        }
        for record in reader {
            match record {
                Ok(record) => {
                    if filter{
                        for name in names.as_ref().unwrap().iter() {
                            if record.name() == name.as_bytes(){
                                let preftag: Vec<u8> = lab_vec[pos].as_bytes().to_vec();
                                let newrecord = edit_record(&record, preftag);
                                match newrecord {
                                    Ok(_) => {
                                        pass_count+=1;
                                        pass_writer.write(&newrecord.unwrap()).unwrap();
                                     },
                                    Err(_) => {
                                        fail_count+=1;
                                        fail_writer.write(&record).unwrap();
                                    },
                                }
                            } else {
                                continue;
                            }
                        }
                    } else {
                        let preftag: Vec<u8> = lab_vec[pos].as_bytes().to_vec();
                        let newrecord  = edit_record(&record, preftag);
                        match newrecord {
                            Ok(_) => {
                                        pass_count+=1;
                                        pass_writer.write(&newrecord.unwrap()).unwrap();
                                     },
                            Err(_) => {
                                fail_count+=1;
                                fail_writer.write(&record).unwrap();
                            },
                        }
                    }
                },
                Err(_) => other_count+=1,
            }
        }
    }
    eprintln!("Processed all reads!!\nFound:\n{} - reads PASSING\n{} - reads PASSING but with issues\n{} - reads FAILING", pass_count, other_count, fail_count);
    return params;
}
    
fn edit_record(record: &bam::record::Record, preftag: Vec<u8>)-> Result<bam::record::Record, Box<dyn Error>>{
    let mut newrecord = record.clone();
    match newrecord.tags().get(b"CB") {
        Some(TagValue::String(array_view, _)) => {
            let oldtag = array_view.to_vec();
            let new_cb = &[preftag, oldtag].concat();
            newrecord.tags_mut().remove(b"CB");
            newrecord.tags_mut().push_string(b"CB", &new_cb);
            Ok(newrecord)
        },
        Some(TagValue::Char(value)) => {
            let oldtag = value.to_string().as_bytes().to_vec();
            let new_cb = &[preftag, oldtag].concat();
            newrecord.tags_mut().remove(b"CB");
            newrecord.tags_mut().push_string(b"CB", &new_cb);
            Ok(newrecord)
        },
        _ => {
            Err(Box::new(std::io::Error::new(std::io::ErrorKind::Other, "'CB' tag not found")))
        }
    }
}

