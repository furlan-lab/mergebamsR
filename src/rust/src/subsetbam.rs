
use std::cmp;
use failure::Error;
use rayon::prelude::*;
use rust_htslib::bam;
use rust_htslib::bam::record::Aux;
use rust_htslib::bam::Record;
use log::{error, info, LevelFilter};
use simplelog::*;
use std::collections::HashMap;
use std::fs;
use std::io::{self, prelude::*};
use std::path::{Path, PathBuf};
use std::process;
use tempfile::tempdir;

pub struct Metrics {
    pub total_reads: usize,
    pub dumped: usize,
    pub kept_reads: usize,
}

pub struct Args<'a> {
    cell_barcodes: &'a HashMap<&'a Vec<u8>, usize>,
    outputbam_no: usize,
    i: usize,
    bam_file: &'a str,
    tmp_dir: &'a Path,
    bam_tag: String,
    virtual_start: Option<i64>,
    virtual_stop: Option<i64>,
    field: &'a str,
    dump_bam: Option<&'a str>,
}


pub struct Outs {
    metrics: Metrics,
    out_paths: Vec<PathBuf>,
}

pub fn subset_bam(
    inputbam: &str,
    final_tags: Vec<Vec<Vec<u8>>>,
    final_outputbams: Vec<String>,
    tag: &str,
    cores: u64,
    field: &str,
    dump_bam: Option<&str>,
) {
    let ll = LevelFilter::Info;
    let bam_tag = tag.to_string();
    let outputbam_no = final_outputbams.len();

    let _ = SimpleLogger::init(ll, Config::default());
    check_inputs_exist(inputbam, final_outputbams.clone());

    let tmp_dir = tempdir().unwrap();
    let virtual_offsets = bgzf_noffsets(inputbam, &cores).unwrap();

    let cell_barcodes: HashMap<_, _> = final_tags
        .iter()
        .enumerate()
        .flat_map(|(index, vec)| vec.iter().map(move |value| (value, index)))
        .collect();

    let chunks: Vec<_> = virtual_offsets
        .iter()
        .enumerate()
        .map(|(i, (virtual_start, virtual_stop))| Args {
            cell_barcodes: &cell_barcodes,
            outputbam_no,
            i,
            bam_file: inputbam,
            tmp_dir: tmp_dir.path(),
            bam_tag: bam_tag.clone(),
            virtual_start: *virtual_start,
            virtual_stop: *virtual_stop,
            field,
            dump_bam,
        })
        .collect();

    let pool = rayon::ThreadPoolBuilder::new()
        .num_threads(cores as usize)
        .build()
        .unwrap();
    let results: Vec<_> = pool.install(|| {
        chunks.par_iter().map(|chunk| slice_bam_chunk(chunk)).collect()
    });

    let mut metrics = Metrics {
        total_reads: 0,
        dumped: 0,
        kept_reads: 0,
    };

    for c in &results {
        add_metrics(&mut metrics, &c.metrics);
    }

    if metrics.kept_reads == 0 {
        error!("Zero alignments were kept. Does your BAM contain the cell barcodes and/or tag you chose?");
        process::exit(1);
    }

    let tmp_bams_vec = transpose_vec(
        results
            .iter()
            .map(|c| &c.out_paths)
            .collect::<Vec<&Vec<PathBuf>>>(),
    );

    if cores == 1 {
        for (i, filefrom) in tmp_bams_vec[0].iter().enumerate() {
            fs::copy(filefrom, &final_outputbams[i]).unwrap();
        }
    } else {
        for (i, tmp_bams) in tmp_bams_vec.into_iter().enumerate() {
            merge_bams(&tmp_bams, Path::new(&final_outputbams[i]));
        }
    }

    info!("Done!");
    info!(
        "Visited {} alignments, dumped {} and kept {}",
        metrics.total_reads, metrics.dumped, metrics.kept_reads
    );
}

fn transpose_vec<T>(v: Vec<&Vec<T>>) -> Vec<Vec<T>>
where
    T: Clone,
{
    assert!(!v.is_empty());
    (0..v[0].len())
        .map(|i| v.iter().map(|inner| inner[i].clone()).collect::<Vec<T>>())
        .collect()
}

fn check_inputs_exist(bam_file: &str, out_bams_path: Vec<String>) {
    if !Path::new(bam_file).exists() {
        error!("File {} does not exist", bam_file);
        process::exit(1);
    }
    for out_bam in out_bams_path.iter() {
        let path = Path::new(out_bam);
        if path.exists() {
            error!("Output path already exists");
            process::exit(1);
        }
        if path.is_dir() {
            error!("Output path is a directory");
            process::exit(1);
        }
        let parent_dir = path.parent().unwrap();
        if !parent_dir.exists() {
            error!("Output directory {:?} does not exist", parent_dir);
            process::exit(1);
        }
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

fn get_tag(rec: &Record, bam_tag: &str) -> Option<Vec<u8>> {
    match rec.aux(bam_tag.as_bytes()) {
        Ok(Aux::String(hp)) => Some(hp.as_bytes().to_vec()),
        _ => None,
    }
}

fn load_writer(bam: &bam::Reader, out_bam_path: &Path) -> Result<bam::Writer, Error> {
    use rust_htslib::bam::Read;
    let hdr = rust_htslib::bam::Header::from_template(bam.header());
    let out_handle = bam::Writer::from_path(out_bam_path, &hdr, rust_htslib::bam::Format::Bam)?;
    Ok(out_handle)
}

fn bgzf_noffsets(
    bam_path: &str,
    num_chunks: &u64,
) -> Result<Vec<(Option<i64>, Option<i64>)>, Error> {
    fn vec_diff(input: &Vec<u64>) -> Vec<u64> {
        let vals = input.iter();
        let next_vals = input.iter().skip(1);
        vals.zip(next_vals).map(|(cur, next)| next - cur).collect()
    }

    if *num_chunks == 1 {
        return Ok(vec![(None, None)]);
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

    let mut final_offsets = Vec::new();
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

fn is_valid_bgzf_block(block: &[u8]) -> bool {
    if block.len() < 18 {
        return false;
    }
    if (block[0] != 31) | (block[1] != 139) | (block[2] != 8) | (block[3] != 4) {
        return false;
    }
    true
}

fn slice_bam_chunk(args: &Args) -> Outs {
    let mut bam = bam::Reader::from_path(args.bam_file).unwrap();
    let mut out_writers: Vec<bam::Writer> = vec![];
    let tmp_out_bam_files: Vec<_> = (0..args.outputbam_no)
        .map(|x| args.tmp_dir.join(format!("tmp_chunk{}_out{}.bam", args.i, x)))
        .collect();

    for tmp_out_bam_file in &tmp_out_bam_files {
        out_writers.push(load_writer(&bam, Path::new(&tmp_out_bam_file)).unwrap());
    }

    let mut dump_writer = if let Some(dump_bam_file) = args.dump_bam {
        Some(load_writer(&bam, Path::new(dump_bam_file)).unwrap())
    } else {
        None
    };

    let mut metrics = Metrics {
        total_reads: 0,
        dumped: 0,
        kept_reads: 0,
    };

    for r in bam.iter_chunk(args.virtual_start, args.virtual_stop) {
        let rec = r.unwrap();
        metrics.total_reads += 1;

        match args.field {
            "name" => {
                let found = process_name(&rec, args, &mut metrics, &mut out_writers);
                if !found & dump_writer.is_some() {
                    metrics.dumped+=1;
                    let _ = dump_writer.as_mut().unwrap().write(&rec);
                }
            }
            "tag" => {
                let found = process_tag(&rec, args, &mut metrics, &mut out_writers);
                if !found & dump_writer.is_some() {
                    metrics.dumped+=1;
                    let _ = dump_writer.as_mut().unwrap().write(&rec);
                }
            }
            _ => {
                error!("Invalid field");
                process::exit(1);
            }
        }
    }

    Outs {
        metrics,
        out_paths: tmp_out_bam_files.clone(),
    }
}

fn process_tag(
    rec: &Record,
    args: &Args,
    metrics: &mut Metrics,
    out_writers: &mut Vec<bam::Writer>,
)  -> bool {
    let barcode = get_tag(&rec, &args.bam_tag).unwrap();
    if let Some(index) = args.cell_barcodes.get(&barcode) {
        metrics.kept_reads += 1;
        out_writers[*index].write(&rec).unwrap();
        return true
    } else {
        return false}
}

fn process_name(
    rec: &Record,
    args: &Args,
    metrics: &mut Metrics,
    out_writers: &mut Vec<bam::Writer>,
)  -> bool {
    let name = rec.qname().to_vec();
    if let Some(index) = args.cell_barcodes.get(&name) {
        metrics.kept_reads += 1;
        out_writers[*index].write(&rec).unwrap();
        return true
    } else {
        return false;
    }
}



fn merge_bams(tmp_bams: &Vec<PathBuf>, out_bam_file: &Path) {
    use bam::Read;
    let bam = bam::Reader::from_path(tmp_bams[0].to_str().unwrap()).unwrap();
    let mut out_bam = load_writer(&bam, out_bam_file).unwrap();
    for b in tmp_bams.iter() {
        let mut rdr = bam::Reader::from_path(b).unwrap();
        for rec in rdr.records() {
            out_bam.write(&rec.unwrap()).unwrap();
        }
    }
}

fn add_metrics(metrics: &mut Metrics, m: &Metrics) {
    metrics.total_reads += m.total_reads;
    metrics.dumped += m.dumped;
    metrics.kept_reads += m.kept_reads;
}




#[cfg(test)]
mod tests {
    use super::*;
    use data_encoding::HEXUPPER;
    use ring::digest::{Context, Digest, SHA256};
    // use tempfile::tempdir;

    pub fn get_library_location()-> String  {
        match std::env::current_exe() {
            Ok(_exe_path) => {
                let path:String = std::env::current_exe().unwrap().to_str().unwrap().to_string();
                let path = path.split("/src/rust/target/debug/deps").collect::<Vec<&str>>()[0];
                let fileloc = Path::new(&path).join("inst/extdata/");
                if fileloc.exists() {
                    fileloc.to_str().unwrap().to_string()
                } else {
                    Path::new(&path).join("extdata/").to_str().unwrap().to_string()
                }
            },
            Err(_) => "failed to get current exe path".to_string()
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

    #[test]
    fn test_bam_single_core() {
        let final_tags = vec![vec![b"ATTGGACAGTCATGCT-1".to_vec(), b"ATCATGGCAGACGCTC-1".to_vec()]];
        let root = get_library_location();
        let out_dir = Path::new(&root).join("test/out");
        fs::create_dir(&out_dir).unwrap();
        let inputbam =  Path::new(&root).join("test/bam1.bam").to_str().unwrap().to_string();
        // let final_outputbams =  Path::new(&root).join("test/out/subset.bam").to_str().unwrap().to_string();
        let final_outputbams1 =  Path::new(&root).join("test/out/subset1_sc.bam").to_str().unwrap().to_string();
        // let final_outputbams2 =  Path::new(&root).join("test/out/subset2.bam").to_str().unwrap().to_string();
        let tag = "CB";
        subset_bam(&inputbam, final_tags, vec![final_outputbams1.clone()], tag, 1, "tag", None);
        let fh = fs::File::open(Path::new(&final_outputbams1)).unwrap();
        let d = sha256_digest(fh).unwrap();
        let d = HEXUPPER.encode(d.as_ref());
        // eprint!("SHA256: {}\n", d);
        assert_eq!(
            d,
            "43F97A078D9860A4D083019919B17CDD05102C71467B223DF201E03B28AB14AD"
        );
        fs::remove_dir_all(out_dir).unwrap();   
    }

    #[test]
    fn test_bam_multiple_core() {
        let final_tags = vec![vec![b"ATTGGACAGTCATGCT-1".to_vec(), b"ATCATGGCAGACGCTC-1".to_vec()],
                                                vec![ b"GGAAAGCTCTCAACTT-1".to_vec(), b"GAGCAGACAGACAGGT-1".to_vec()]];
        let root = get_library_location();
        let out_dir = Path::new(&root).join("test/out");
        fs::create_dir(&out_dir).unwrap();
        let inputbam =  Path::new(&root).join("test/bam1.bam").to_str().unwrap().to_string();
        let final_outputbams1 =  Path::new(&root).join("test/out/subset1.bam").to_str().unwrap().to_string();
        let final_outputbams2 =  Path::new(&root).join("test/out/subset2.bam").to_str().unwrap().to_string();
        let tag = "CB";
        subset_bam(&inputbam, final_tags, vec![final_outputbams1.clone(), final_outputbams2.clone()], tag, 8, "tag", None);
        let fh = fs::File::open(Path::new(&final_outputbams2)).unwrap();
        let d = sha256_digest(fh).unwrap();
        let d = HEXUPPER.encode(d.as_ref());
        // eprint!("SHA256: {}\n", d);
        assert_eq!(
            d,
            "536FC3AEEF2B007AAAB85B7E29F8E6BDCBEAF7A8F1D14C9A0EC849FABCDC1E1D"
        );
        let fh = fs::File::open(Path::new(&final_outputbams1)).unwrap();
        let d = sha256_digest(fh).unwrap();
        let d = HEXUPPER.encode(d.as_ref());
        assert_eq!(
            d,
            "43F97A078D9860A4D083019919B17CDD05102C71467B223DF201E03B28AB14AD"
        );
        fs::remove_dir_all(out_dir).unwrap();  
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