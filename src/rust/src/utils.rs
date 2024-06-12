extern crate bam;
extern crate csv;
extern crate flate2;
extern crate itertools;

// use bam::record;
// use itertools::Itertools;
// use differ::{Differ, Tag};
// use bam::RecordWriter;
use bam::record::tags::TagValue;
use std::str;
use std::error::Error;

// pub fn subset_bam_rust(inputbam: &str, final_tags: Vec<Vec<String>>, outputbams: Vec<String>, prefixes: Vec<String>, tag: &str) {
//     let reader = bam::BamReader::from_path(inputbam.to_string(), 0).unwrap();
//     let mut writers: Vec<bam::BamWriter<Writer>> = Vec::new();
//     for outputbam in outputbams {
//         let writer = bam::BamWriter::from_path(outputbam.to_string(), reader.header().clone()).unwrap();
//         writers.push(writer);
//     }
//     let mut pass_count: u64 = 0;
//     let mut fail_count: u64 = 0;
//     for (record_num, record) in reader.enumerate() {
//         match record {
//             Ok(record) => {
//                 let tag_value = get_tag(record, tag);
//                 match tag_value {
//                     Ok(tag_value) => {
//                         let tag_value_str = str::from_utf8(&tag_value).unwrap();
//                         for (i, final_tag) in final_tags.iter().enumerate() {
//                             if final_tag.contains(&tag_value_str.to_string()) {
//                                 writers[i].write(&record).unwrap();
//                             }
//                         }
//                         pass_count+=1;
//                     },
//                     Err(_) => fail_count+=1,
//                 }
//             },
//             Err(_) => fail_count+=1,
//         }
//     }
//     eprint!("Read {} records\n", pass_count);
// }

pub fn peekbam_rust<'a>(bam: &'a str, n: u64, field: &'a str, tag: &'a str) -> Result<Vec<String>, extendr_api::Error> {
    let reader = bam::BamReader::from_path(bam.to_string(), 0).unwrap();
    let mut pass_count: u64 = 0;
    let mut fail_count: u64 = 0;
    let mut tags: Vec<String> = Vec::new();
    match field {
        "tag" => {
            for record in reader {
                match record {
                    Ok(record) => {
                        let tag_value = get_tag(record, tag);
                        match tag_value {
                            Ok(tag_value) => {
                                if pass_count < n {
                                    tags.push(String::from_utf8(tag_value).unwrap());
                                    pass_count+=1;
                                } else {
                                    break;
                                }
                            },
                            Err(_) => fail_count+=1,
                        }
                    },
                    Err(_) => fail_count+=1,
                }
            }
        },
        "name" => {
            for record in reader {
                match record {
                    Ok(record) => {
                        let value: Result<Vec<u8>, extendr_api::Error> = Ok(record.name().to_vec());
                        match value {
                            Ok(tag_value) => {
                                if pass_count < n {
                                    tags.push(String::from_utf8(tag_value).unwrap());
                                    pass_count+=1;
                                } else {
                                    break;
                                }
                            },
                            Err(_) => fail_count+=1,
                        }
                    },
                    Err(_) => fail_count+=1,
                }
            }
        },
        _ => {
            eprint!("ERROR: field not recognized\n");
            return Err(extendr_api::Error::from("field not recognized"));
        }
    }
    
    eprint!("Read {} records\n", pass_count);
    eprint!("Failed to read {} records\n", fail_count);
    Ok(tags)
}

fn get_tag(record: bam::Record, tag: &str) -> Result<Vec<u8>, Box<dyn Error>> {
    let tag_value = record.tags().get(&[tag.as_bytes()[0], tag.as_bytes()[1]]);
    match tag_value {
        Some(TagValue::String(array_view, _)) => {
            let captured_tag = array_view.to_vec();
            Ok(captured_tag)
        },
        Some(TagValue::Char(value)) => {
            let captured_tag = value.to_string().as_bytes().to_vec();
            Ok(captured_tag)
        },        
        _ => {
            Err(Box::new(std::io::Error::new(std::io::ErrorKind::Other, "Tag not found")))
        }
    }
}
