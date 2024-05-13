#![allow(unused_variables)]
// extern crate bam;
// use clap::{App, Arg};
// use faccess::{AccessMode, PathExt};
use failure::Error;
use rayon::prelude::*;
use rust_htslib::bam;
use rust_htslib::bam::record::Aux;
use rust_htslib::bam::Record;
// use rust_htslib::bam::Read;
// use bam::record::{Aux, Record};
// use bam::record::Aux;
// use bam::Record;
// use bam::BamReader;
// use bam::BamWriter;
use log::error;
use log::info;
use simplelog::*;
use std::cmp;
use std::collections::HashSet;
use std::fs;
use std::io::prelude::*;
use std::io::{self};
use std::path::{Path, PathBuf};
use std::process;
use tempfile::tempdir;
// use std::collections::HashMap;
use std::collections::HashMap;

pub struct Metrics {
    pub total_reads: usize,
    pub barcoded_reads: usize,
    pub kept_reads: usize,
}

// pub struct ChunkArgs<'a> {
//     cell_barcodes: &'a HashSet<Vec<u8>>,
//     i: usize,
//     bam_file: &'a str,
//     tmp_dir: &'a Path,
//     bam_tag: String,
//     virtual_start: Option<i64>,
//     virtual_stop: Option<i64>,
// }
pub struct ChunkArgs<'a> {
    cell_barcodes: &'a HashMap<&'a Vec<u8>, usize>,
    i: usize,
    bam_file: &'a str,
    // bam_file: Vec<String>,
    tmp_dir: Vec<String>,
    bam_tag: String,
    virtual_start: Option<i64>,
    virtual_stop: Option<i64>,
}


pub struct BamArgs<'a> {
    cell_barcodes: &'a HashSet<Vec<u8>>,
    bam_file: &'a str,
    out_dir: &'a Path,
    bam_tag: String,
}

// pub struct ChunkOuts {
//     metrics: Metrics,
//     out_bam_file: PathBuf,
// }

pub struct ChunkOuts {
    metrics: Metrics,
    out_paths: Vec<PathBuf>,
}


// pub fn subset_bam_rust(inputbam: &str, final_tags: Vec<Vec<u8>>, final_outputbams: Vec<String>, final_prefixes: Vec<String>, tag: &str) {
//     let ll = "error";
//     let bam_tag = tag.to_string();
//     let out_bam_file = &final_outputbams[0].to_string();
//     let ll = match ll {
//         "info" => LevelFilter::Info,
//         "debug" => LevelFilter::Debug,
//         "error" => LevelFilter::Error,
//         &_ => {
//             eprintln!("Log level not valid");
//             process::exit(1);
//         }
//     };
//     let _ = SimpleLogger::init(ll, Config::default());
//     let bam_file = inputbam;
    
//     check_inputs_exist(bam_file, out_bam_file);
//     let cell_barcodes = final_tags.iter().cloned().collect();
//     let c = BamArgs {
//         cell_barcodes: &cell_barcodes,
//         bam_file: &bam_file,
//         out_dir: Path::new(out_bam_file),
//         bam_tag: bam_tag.clone(),
//     };

//     let results = slice_bam(&c);

//     // combine metrics
//     let mut metrics = Metrics {
//         total_reads: 0,
//         barcoded_reads: 0,
//         kept_reads: 0,
//     };

//     fn add_metrics(metrics: &mut Metrics, m: &Metrics) {
//         metrics.total_reads += m.total_reads;
//         metrics.barcoded_reads += m.barcoded_reads;
//         metrics.kept_reads += m.kept_reads;
//     }

//     // let mut tmp_bams = Vec::new();
//     add_metrics(&mut metrics, &results.metrics);

//     if metrics.kept_reads == 0 {
//         error!("Zero alignments were kept. Does your BAM contain the cell barcodes and/or tag you chose?");
//         process::exit(2);
//     }

//     info!("Done!");
//     info!(
//         "Visited {} alignments, found {} with barcodes and kept {}",
//         metrics.total_reads, metrics.barcoded_reads, metrics.kept_reads
//     );
// }



pub fn subset_bam_rust_parallel(inputbam: &str, final_tags: Vec<Vec<u8>>, final_outputbams: Vec<String>, final_prefixes: Vec<String>, tag: &str, cores: u64) {
    
    // use std::collections::HashMap;
    let ll = "error";
    let bam_tag = tag.to_string();
    // let out_bam_file = &final_outputbams[0].to_string();
    let ll = match ll {
        "info" => LevelFilter::Info,
        "debug" => LevelFilter::Debug,
        "error" => LevelFilter::Error,
        &_ => {
            eprintln!("Log level not valid");
            process::exit(1);
        }
    };
    let _ = SimpleLogger::init(ll, Config::default());
    let bam_file = inputbam;
    
    // check_inputs_exist(bam_file, out_bam_file);
    // let cell_barcodes = final_tags.iter().cloned().collect();
    let tmp_dir = tempdir().unwrap();
    let virtual_offsets = bgzf_noffsets(&bam_file, &cores).unwrap();
    let mut cell_barcodes = HashMap::new();
    for vec in final_tags.iter().enumerate(){
        cell_barcodes.insert(vec.1, vec.0);
        // for cb in vec.1.iter() {
        //     cell_barcodes.insert(cb, vec.0);
        // }
    }

    let mut chunks = Vec::new();
    for (i, (virtual_start, virtual_stop)) in virtual_offsets.iter().enumerate() {
        let c = ChunkArgs {
            cell_barcodes: &cell_barcodes,
            i: i,
            bam_file: bam_file,
            tmp_dir: final_outputbams.clone(),
            bam_tag: bam_tag.clone(),
            virtual_start: *virtual_start,
            virtual_stop: *virtual_stop,
        };
        chunks.push(c);
    }
    let pool = rayon::ThreadPoolBuilder::new()
        .num_threads(cores as usize)
        .build()
        .unwrap();
    let results: Vec<_> = pool.install(|| {
        chunks
            .par_iter()
            .map(|chunk| slice_bam_chunk(chunk))
            .collect()
    });

    // combine metrics
    let mut metrics = Metrics {
        total_reads: 0,
        barcoded_reads: 0,
        kept_reads: 0,
    };

    fn add_metrics(metrics: &mut Metrics, m: &Metrics) {
        metrics.total_reads += m.total_reads;
        metrics.barcoded_reads += m.barcoded_reads;
        metrics.kept_reads += m.kept_reads;
    }

    let mut tmp_bams = Vec::new();
    for c in results.iter() {
            add_metrics(&mut metrics, &c.metrics);
            tmp_bams.push(&c.out_paths);
    }

    // if metrics.kept_reads == 0 {
    //     error!("Zero alignments were kept. Does your BAM contain the cell barcodes and/or tag you chose?");
    //     process::exit(2);
    // }

    // just copy the temp file over
    for out_bam_file in final_outputbams{
        if cores == 1 {
            fs::copy(tmp_bams[0], out_bam_file).unwrap();
        } else {
            info!("Merging {} BAM chunks into final output", cores);
            merge_bams(tmp_bams.clone(), Path::new(&out_bam_file));
        }
    }


    info!("Done!");
    info!(
        "Visited {} alignments, found {} with barcodes and kept {}",
        metrics.total_reads, metrics.barcoded_reads, metrics.kept_reads
    );
}

pub fn check_inputs_exist(bam_file: &str, out_bam_path: &str) {
    if !Path::new(&bam_file).exists() {
        error!("File {} does not exist", bam_file);
        process::exit(1);
    }
    let path = Path::new(out_bam_path);
    if path.exists() {
        error!("Output path already exists");
        process::exit(1);
    }
    if path.is_dir() {
        error!("Output path is a directory");
        process::exit(1);
    }
    let _parent_dir = path.parent();
    if _parent_dir.is_none() {
        error!("Unable to parse directory from {}", out_bam_path);
        process::exit(1);
    }
    let parent_dir = _parent_dir.unwrap();
    if (parent_dir.to_str().unwrap().len() > 0) & !parent_dir.exists() {
        error!("Output directory {:?} does not exist", parent_dir);
        process::exit(1);
    }

    let extension = Path::new(bam_file).extension().unwrap().to_str().unwrap();
    match extension {
        "bam" => {
            let bai = bam_file.to_owned() + ".bai";
            if !Path::new(&bai).exists() {
                error!("BAM index {} does not exist", bai);
                process::exit(1);
            }
        }
        "cram" => {
            let crai = bam_file.to_owned() + ".crai";
            if !Path::new(&crai).exists() {
                error!("CRAM index {} does not exist", crai);
                process::exit(1);
            }
        }
        &_ => {
            error!("BAM file did not end in .bam or .cram. Unable to validate");
            process::exit(1);
        }
    }
}

pub fn get_cell_barcode(rec: &Record, bam_tag: &str) -> Option<Vec<u8>> {
    // println!("{:?}", rec.aux(bam_tag.as_bytes()));
    match rec.aux(bam_tag.as_bytes()) {
        Ok(Aux::String(hp)) => {
            let cb = hp.as_bytes().to_vec();
            Some(cb)
        }
        _ => None,
    }
}

// pub fn load_writer(bam: &bam::BamReader<, out_bam_path: &Path) -> Result<bam::RecordWriter, Error> {
//     use rust_htslib::bam::Read; // collides with fs::Read
//     let hdr = rust_htslib::bam::Header::from_template(bam.header());
//     let out_handle = bam::RecordWriter::from_path(out_bam_path, &hdr, rust_htslib::bam::Format::Bam)?;
//     Ok(out_handle)
// }

pub fn load_writer(bam: &bam::Reader, out_bam_path: &Path) -> Result<bam::Writer, Error> {
    use rust_htslib::bam::Read; // collides with fs::Read
    let hdr = rust_htslib::bam::Header::from_template(bam.header());
    let out_handle = bam::Writer::from_path(out_bam_path, &hdr, rust_htslib::bam::Format::Bam)?;
    Ok(out_handle)
}

pub fn bgzf_noffsets(
    bam_path: &str,
    num_chunks: &u64,
) -> Result<Vec<(Option<i64>, Option<i64>)>, Error> {
    fn vec_diff(input: &Vec<u64>) -> Vec<u64> {
        let vals = input.iter();
        let next_vals = input.iter().skip(1);

        vals.zip(next_vals).map(|(cur, next)| next - cur).collect()
    }

    // if we only have one thread, this is easy
    if *num_chunks == 1 as u64 {
        let final_offsets = vec![(None, None)];
        return Ok(final_offsets);
    }

    let bam_bytes = fs::metadata(bam_path)?.len();
    let mut initial_offsets = Vec::new();
    let step_size = bam_bytes / num_chunks;
    for n in 1..*num_chunks {
        initial_offsets.push((step_size * n) as u64);
    }

    let num_bytes = if initial_offsets.len() > 1 {
        let diff = vec_diff(&initial_offsets);
        let m = diff.iter().max().unwrap();
        cmp::min(1 << 16, *m)
    } else {
        1 << 16
    };

    // linear search to the right of each possible offset until
    // a valid virtual offset is found
    let mut adjusted_offsets = Vec::new();
    let mut fp = fs::File::open(bam_path)?;
    for offset in initial_offsets {
        fp.seek(io::SeekFrom::Start(offset))?;
        let mut buffer = [0; 2 << 16];
        fp.read(&mut buffer)?;
        for i in 0..num_bytes {
            if is_valid_bgzf_block(&buffer[i as usize..]) {
                adjusted_offsets.push(offset + i);
                break;
            }
        }
    }
    // bit-shift and produce start/stop intervals
    let mut final_offsets = Vec::new();

    // handle special case where we only found one offset
    if adjusted_offsets.len() == 1 {
        final_offsets.push((None, None));
        return Ok(final_offsets);
    }

    final_offsets.push((None, Some(((adjusted_offsets[1]) as i64) << 16)));
    for n in 2..num_chunks - 1 {
        let n = n as usize;
        final_offsets.push((
            Some((adjusted_offsets[n - 1] as i64) << 16),
            Some((adjusted_offsets[n] as i64) << 16),
        ));
    }
    final_offsets.push((
        Some(((adjusted_offsets[adjusted_offsets.len() - 1]) as i64) << 16),
        None,
    ));
    Ok(final_offsets)
}

pub fn is_valid_bgzf_block(block: &[u8]) -> bool {
    // look for the bgzip magic characters \x1f\x8b\x08\x04
    // TODO: is this sufficient?
    if block.len() < 18 {
        return false;
    }
    if (block[0] != 31) | (block[1] != 139) | (block[2] != 8) | (block[3] != 4) {
        return false;
    }
    true
}

pub fn slice_bam_chunk(args: &ChunkArgs) -> ChunkOuts {
    let mut bam = bam::Reader::from_path(args.bam_file).unwrap();
    // let out_bam_file = args.tmp_dir.join(format!("{}.bam", args.i));
    let out_bam_files = &args.tmp_dir;
    let mut out_writers: Vec<bam::Writer> =  vec![];
    for out_bam_file in out_bam_files{
        out_writers.push(load_writer(&bam, Path::new(&out_bam_file)).unwrap());
    }
    // let mut out_bam = load_writer(&bam, &out_bam_file).unwrap();
    let mut metrics = Metrics {
        total_reads: 0,
        barcoded_reads: 0,
        kept_reads: 0,
    };
    for r in bam.iter_chunk(args.virtual_start, args.virtual_stop) {
        let rec = r.unwrap();
        metrics.total_reads += 1;
        let barcode = get_cell_barcode(&rec, &args.bam_tag);
        if barcode.is_some() {
            metrics.barcoded_reads += 1;
            let barcode = barcode.unwrap();
            let out_bam_index = args.cell_barcodes.get(&barcode);
            match out_bam_index{
                Some(index) => {metrics.kept_reads +=1;
                                        out_writers[*index].write(&rec).unwrap();
                                    }
                None => continue
            }
            // if args.cell_barcodes.contains(&barcode) {
            //     metrics.kept_reads += 1;
            //     out_bam.write(&rec).unwrap();
            // }
        }
    }
    let mut out_paths: Vec<_> = vec![];

    for out_bam_file in out_bam_files{
        out_paths.push();
    }
    let r = ChunkOuts {
        metrics: metrics,
        // out_bam_file: Path::new(&out_bam_file).to_path_buf(),
        // out_bam_file: Path::new(&out_bam_file).to_path_buf(),
        out_paths: out_paths,
    };
    info!("Chunk {} is done", args.i);
    r
}



// pub fn slice_bam(args: &BamArgs) -> ChunkOuts {
//     use crate::rust_htslib::bam::Read;
//     let mut bam = bam::Reader::from_path(args.bam_file).unwrap();
//     // let out_bam_file = args.out_dir.join(format!("{}.bam", args.i));
//     let mut out_bam = load_writer(&bam, &args.out_dir).unwrap();
//     let mut metrics = Metrics {
//         total_reads: 0,
//         barcoded_reads: 0,
//         kept_reads: 0,
//     };
//     for r in bam.records(){
//         let rec = r.unwrap();
//         metrics.total_reads += 1;
//         let barcode = get_cell_barcode(&rec, &args.bam_tag);
//         if barcode.is_some() {
//             metrics.barcoded_reads += 1;
//             let barcode = barcode.unwrap();
//             if args.cell_barcodes.contains(&barcode) {
//                 metrics.kept_reads += 1;
//                 out_bam.write(&rec).unwrap();
//             }
//         }
//     }
//     let r = ChunkOuts {
//         metrics: metrics,
//         out_bam_file: args.out_dir.to_path_buf(),
//     };
//     r
// }



pub fn merge_bams(tmp_bams: Vec<&PathBuf>, out_bam_file: &Path) {
    use rust_htslib::bam::Read; // collides with fs::Read
    let bam = bam::Reader::from_path(tmp_bams[0]).unwrap();
    let mut out_bam = load_writer(&bam, out_bam_file).unwrap();
    for b in tmp_bams.iter() {
        let mut rdr = bam::Reader::from_path(b).unwrap();
        for _rec in rdr.records() {
            let rec = _rec.unwrap();
            out_bam.write(&rec).unwrap();
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use data_encoding::HEXUPPER;
    use ring::digest::{Context, Digest, SHA256};
    // use tempfile::tempdir;

    fn get_library_location()-> String  {
        match std::env::current_exe() {
            Ok(exe_path) => {
                let path:String = std::env::current_exe().unwrap().to_str().unwrap().to_string();
                let path = path.split("/src/rust/target/debug/deps").collect::<Vec<&str>>()[0];
                let fileloc = Path::new(&path).join("inst/extdata/");
                if fileloc.exists() {
                    fileloc.to_str().unwrap().to_string()
                } else {
                    Path::new(&path).join("extdata/").to_str().unwrap().to_string()
                }
            },
            Err(e) => "failed to get current exe path".to_string()
        }
    }

    /// Compute digest value for given `Reader` and print it
    /// This is taken from the Rust cookbook
    /// https://rust-lang-nursery.github.io/rust-cookbook/cryptography/hashing.html
    fn sha256_digest<R: Read>(mut reader: R) -> Result<Digest, Error> {
        let mut context = Context::new(&SHA256);
        let mut buffer = [0; 1024];

        loop {
            let count = reader.read(&mut buffer)?;
            if count == 0 {
                break;
            }
            context.update(&buffer[..count]);
        }

        Ok(context.finish())
    }

    // #[test]
    // fn test_bam_single_core() {
    //     let final_tags = vec![b"ATTGGACAGTCATGCT-1".to_vec(), b"TTTACTGAGTCGATAA-1".to_vec()];
    //     let root = get_library_location();
    //     let inputbam =  Path::new(&root).join("test/bam1.bam").to_str().unwrap().to_string();
    //     let final_prefixes = vec!["".to_string()];
    //     let final_outputbams =  Path::new(&root).join("test/out/subset.bam").to_str().unwrap().to_string();
    //     let tag = "CB";
    //     subset_bam_rust(&inputbam, final_tags, vec![final_outputbams.clone()], final_prefixes, tag);
    //     let fh = fs::File::open(Path::new(&final_outputbams)).unwrap();
    //     let d = sha256_digest(fh).unwrap();
    //     let d = HEXUPPER.encode(d.as_ref());
    //     // eprint!("SHA256: {}\n", d);
    //     assert_eq!(
    //         d,
    //         "9D5A7360049F4BB840D6FE2EAE0FE87B708C0816FD7F7F27A919F404D0C966E1"
    //     );
    // }

    #[test]
    fn test_bam_multiple_core() {
        let final_tags = vec![b"ATTGGACAGTCATGCT-1".to_vec(), b"TTTACTGAGTCGATAA-1".to_vec()];
        let root = get_library_location();
        let inputbam =  Path::new(&root).join("test/bam1.bam").to_str().unwrap().to_string();
        let final_prefixes = vec!["".to_string()];
        let final_outputbams =  Path::new(&root).join("test/out/subset.bam").to_str().unwrap().to_string();
        let tag = "CB";
        subset_bam_rust_parallel(&inputbam, final_tags, vec![final_outputbams.clone()], final_prefixes, tag, 4);
        let fh = fs::File::open(Path::new(&final_outputbams)).unwrap();
        let d = sha256_digest(fh).unwrap();
        let d = HEXUPPER.encode(d.as_ref());
        // eprint!("SHA256: {}\n", d);
        assert_eq!(
            d,
            "9D5A7360049F4BB840D6FE2EAE0FE87B708C0816FD7F7F27A919F404D0C966E1"
        );
    }

    #[test]
    fn test_hashmaps() {
    use std::collections::HashMap;
    let mut cells = HashMap::new();

    // Review some books.
    cells.insert(
        "ATTGGACAGTCATGCT-1".to_string(),
        "1".to_string(),
    );
    cells.insert(
        "TTTACTGAGTCGATAA-1".to_string(),
        "1".to_string(),
    );
    cells.insert(
        "ATCATGGCAGACGCTC-1".to_string(),
        "2".to_string(),
    );
    cells.insert(
        "TCATTACTCGGAAACG-1".to_string(),
        "2".to_string(),
    );
    cells.insert(
        "ATGTGTGCACATGTGT-1".to_string(),
        "2".to_string(),
    );
    cells.insert(
        "TTTACTGAGTCGATAA-1".to_string(),
        "3".to_string(),
    );
    cells.insert(
        "AAGTCTGCACAGGAGT-1".to_string(),
        "3".to_string(),
    );
    // eprintln!("{:?}", cells);
    // let ans = cells.get("AAGTCTGCACAGGAGT-1");
    // assert_eq!(ans.unwrap(), "3");
    let mut cells = HashMap::new();
    let tags = vec![vec!["ATTGGACAGTCATGCT-1", "TTTACTGAGTCGATAA-1", "ATCATGGCAGACGCTC-1", "CACTCCATCTCGCTTG-1", "GCTGCGAGTCCGTCAG-1", "CTCTACGGTCCAGTAT-1", "TAGTTGGGTTCAGACT-1", "CAGTCCTAGCAATATG-1"],
                                    vec![ "CGACCTTTCCACGTGG-1", "GGAAAGCTCTCAACTT-1", "GGACATTCAGCAGTTT-1", "GACGCGTCATCTCGCT-1", "TCGGGACGTGCTCTTC-1", "GTTCATTTCGAATCCA-1", "GTGCAGCGTACACCGC-1", "CAGTCCTTCACGGTTA-1"],
                                    vec!["GAGCAGACAGACAGGT-1", "ATGTGTGCACATGTGT-1", "TCATTACTCGGAAACG-1", "TGCGTGGAGCCATCGC-1", "GTCTCGTCACCTATCC-1", "CATTATCCATTGTGCA-1", "AAGTCTGCACAGGAGT-1", "GACCAATGTCCTCTTG-1"]];
    for vec in tags.iter().enumerate(){
        for cb in vec.1.iter() {
            cells.insert(cb, vec.0);
        }
    }

    eprintln!("{:?}", cells.get(&"GAGCAGACAGACAGGT-1").unwrap());
    }
}
