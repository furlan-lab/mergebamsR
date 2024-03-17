#![allow(non_snake_case)]

use extendr_api::prelude::*;
mod mergebamsR;


/// mergebams_rust
/// @export
/// @keywords internal
#[extendr]
// fn mergebams_rust_helper(bams: Robj, out_path: Robj, names: Robj, prefixes: Robj) {
//   let bam_files: Vec<&str> = match bams.as_str_vector() {
//     Some(files) => files,
//     None => return,
//   };
//   let out_path: &str = match out_path.as_str_vector() {
//     Some(files) => files[0],
//     None => return,
//   };
//   let mut read_names: Vec<Vec<&str>> = vec![];
//   names.as_list().unwrap().iter().for_each(|x | {
//     if x.1.is_null(){
//       read_names.push(vec![&""]);
//     } else {
//       let y = x.1.as_str_vector().unwrap();
//       read_names.push(y);
//     }
//     // read_names.push(x.1.as_str_vector().unwrap());
//   });
//   let prefixes: Vec<&str> = prefixes.as_str_vector().unwrap();
//   mergebamsR::mergebams_rust(bam_files, out_path, read_names, prefixes)
// }

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
    mergebamsR::mergebams_rust(bam_files, &out_path, read_names, prefixes);
}


// Macro to generate exports.
// This ensures exported functions are registered with R.
// See corresponding C code in `entrypoint.c`.
extendr_module! {
  mod mergebamsR;
  fn mergebams_rust_helper;
}
