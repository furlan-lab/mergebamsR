extern crate bam;
extern crate csv;
extern crate flate2;
extern crate itertools;

// use itertools::Itertools;
// use differ::{Differ, Tag};
// use bam::RecordWriter;
use bam::record::tags::TagValue;
use std::str;
use std::error::Error;

pub fn peekbam_rust<'a>(bam: &'a str, n: u64, tag: &'a str) -> Result<Vec<String>, extendr_api::Error> {
    let reader = bam::BamReader::from_path(bam.to_string(), 0).unwrap();
    let mut pass_count: u64 = 0;
    let mut fail_count: u64 = 0;
    let mut tags: Vec<String> = Vec::new();
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
    eprint!("Read {} records", pass_count);
    eprint!("Failed to read {} records", fail_count);
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
