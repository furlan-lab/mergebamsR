#![allow(non_snake_case)]

// use std::usize;

// use bam::record::tags;
use extendr_api::prelude::*;
mod mergebams;
mod utils;


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
fn peekbam_rust_helper(bam: Robj, n: Robj, tag: Robj) -> Robj{
    let bam_file: &str  = match bam.as_str_vector() {
        Some(files) => files[0],
        None => {
                  eprintln!("ERROR: bams is not a string vector");
                  return Robj::from(0)
                },
    };
    let n:u64 = match n.as_integer() {
        Some(n) => <i32 as extendr_api::TryInto<u64>>::try_into(n).unwrap(),
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
    let tags = utils::peekbam_rust(bam_file, n, tag);
    // match tags {
    //     Ok(tags) => Robj::from(tags),
    //     Err(_) => {
    //         eprintln!("ERROR: peekbam_rust failed");
    //         return Robj::from(0)
    //     }
    // }
    Robj::from(&tags.unwrap())
}


// Macro to generate exports.
// This ensures exported functions are registered with R.
// See corresponding C code in `entrypoint.c`.
extendr_module! {
    mod mergebamsR;
    fn mergebams_rust_helper;
    fn peekbam_rust_helper;
}
