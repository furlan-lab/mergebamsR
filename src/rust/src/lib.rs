#![allow(non_snake_case)]
#![allow(unused_assignments)]


extern crate clap;
extern crate csv;
extern crate data_encoding;
extern crate failure;
extern crate rayon;
extern crate ring;
extern crate rust_htslib;
extern crate simplelog;
extern crate tempfile;
extern crate terminal_size;
extern crate log;
extern crate faccess;
extern crate human_panic;

// use std::usize;

// use bam::record::tags;
use extendr_api::prelude::*;
// use rust_htslib::htslib::sam_fields;
mod mergebams;
mod utils;
mod subsetbam;


/// mergebams_rust
/// @export
/// @keywords internal
#[extendr]
fn mergebams_rust_helper(bams: Robj, out_path: Robj, names: Robj, prefixes: Robj) {
    let bam_files: Vec<&str> = match bams.as_str_vector() {
        Some(files) => files,
        None => {
                  eprintln!("ERROR: bams is not a string vector");
                  return
                },
    };
    let out_path: &str = match out_path.as_str_vector() {
        Some(paths) => paths[0],
        None => {
                  eprintln!("ERROR: out_path is not a string");
                  return
                },
    };
    let prefixes: Vec<&str> = match prefixes.as_str_vector() {
        Some(prefs) => prefs,
        None => {
                  eprintln!("ERROR: prefixes is not a string vector");
                  return
                },
    };
    let mut read_names: Vec<Option<Vec<String>>> = Vec::new();
    if let Some(names_list) = names.as_list() {
        for name in names_list.iter() {
            if name.1.is_null() {
                read_names.push(None);
            } else {
                let names_vec: Vec<String> = name.1.as_str_vector().unwrap_or_default().into_iter().map(String::from).collect();
                read_names.push(Some(names_vec));
            }
        }
    }
    // let prefixes: Vec<&str> = prefixes.as_str_vector().unwrap();

    // Assuming mergebamsR::mergebams_rust now accepts Vec<String> instead of Vec<&str>
    mergebams::mergebams_rust(bam_files, &out_path, read_names, prefixes);
}

/// peekbam_rust
/// @export
/// @keywords internal
#[extendr]
fn peekbam_rust_helper(bam: Robj, n: Robj, field: Robj, tag: Robj) -> Robj{
    let field: &str = match field.as_str_vector() {
        Some(fields) => fields[0],
        None => {
                  eprintln!("ERROR: field is not a string");
                  return Robj::from(0)
                },
    };

    let bam_file: &str  = match bam.as_str_vector() {
        Some(files) => files[0],
        None => {
                  eprintln!("ERROR: bams is not a string vector");
                  return Robj::from(0)
                },
    };
    let n = match n.as_real() {
        Some(n) => n as u64,
        None => {
                  eprintln!("ERROR: n is not an integer");
                  return Robj::from(0)
                },
        };


    let tag: &str = match tag.as_str_vector() {
        Some(tags) => tags[0],
        None => {
                  eprintln!("ERROR: tag is not a string");
                  return Robj::from(0)
                },
    };

    // let mut tags: Result<Vec<&str>> = Ok(Vec::new());
    let tags = utils::peekbam_rust(bam_file, n, field, tag);
    match tags {
        Ok(tags) => Robj::from(tags),
        Err(_) => {
            eprintln!("ERROR: peekbam_rust failed");
            return Robj::from(0)
        }
    }
    // Robj::from(&tags.unwrap())
}

/// subsetbam_rust
/// @export
/// @keywords internal
#[extendr]
fn subsetbam_rust_helper(inputbam: Robj, features: Robj, outputbams: Robj, tag: Robj, cores: Robj, field: Robj, dump_bam: Robj){
    
    let inputbam: &str  = match inputbam.as_str_vector() {
        Some(files) => files[0],
        None => {
                  eprintln!("ERROR: inputbam is not a string");
                  return
                },
    };

    let tag: &str = match tag.as_str_vector() {
        Some(tags) => tags[0],
        None => {
                  eprintln!("ERROR: tag is not a string");
                  return 
                },
    };

    let field: &str = match field.as_str_vector() {
        Some(fields) => fields[0],
        None => {
                  eprintln!("ERROR: field is not a string");
                  return 
                },
    };

    let mut final_features: Vec<Vec<Vec<u8>>> = Vec::new();
    let features_unlisted = features.as_list().unwrap();
    // eprintln!("{:?} tags unlisted length", tags_unlisted.len());
    for (_item_str, item_robj) in features_unlisted {
        let mut int_vec: Vec<Vec<u8>> = Vec::new();
        let data = item_robj.as_string_vector().unwrap();
        for element in data {
            let datain: Vec<u8> = element.as_bytes().to_vec();
            int_vec.push(datain);
        }
        final_features.push(int_vec);
    }
    let final_outputbams = match outputbams.as_string_vector() {
        Some(files) => files,
        None => {
                  eprintln!("ERROR: outputbams is not a string vector");
                  return
                },
    };
    let cores = match cores.as_real() {
        Some(n) => n as u64,
        None => {
                  eprintln!("ERROR: cores is not an integer");
                  return
                },
        };
    let mut dump_bam_r: Option<&str> = None;

    if Robj::is_na(&dump_bam){
        dump_bam_r = None;
    } else {
        dump_bam_r = match dump_bam.as_str_vector() {
            Some(dump) => Some(dump[0]),
            None => {
                      eprintln!("ERROR: dump_bam is not a string");
                      return
                    },
        };
    }

    // if cores>1{
    //     subsetbam::subset_bam_rust_split(inputbam, final_features, final_outputbams, final_prefixes, tag, cores, field, dump_bam_r);
    // } else {
    //     subsetbam::subset_bam_rust(inputbam, final_features, final_outputbams, final_prefixes, tag, field, dump_bam_r);
    // }
    subsetbam::subset_bam(inputbam, final_features, final_outputbams, tag, cores, field, dump_bam_r);
    // subsetbam::subset_bam_rust_parallel(inputbam, final_tags, final_outputbams, final_prefixes, tag, cores);
    
}



// Macro to generate exports.
// This ensures exported functions are registered with R.
// See corresponding C code in `entrypoint.c`.
extendr_module! {
    mod mergebamsR;
    fn mergebams_rust_helper;
    fn peekbam_rust_helper;
    fn subsetbam_rust_helper;
}
