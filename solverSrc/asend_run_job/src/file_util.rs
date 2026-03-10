use crate::cpp_str::CppStr;
use crate::constants::*;

use std::fs::File;
use std::io::{self, Read, Write, BufRead};
use std::path::Path;

pub fn read_lines<P>(filename: P) -> io::Result<io::Lines<io::BufReader<File>>> where P: AsRef<Path>, {
    let file = File::open(filename)?;
    Ok(io::BufReader::new(file).lines())
}

pub fn read_input_line(file_line : &mut CppStr, headings : &mut Vec<CppStr>, hd_ld_space : &mut [usize], data : &mut Vec<CppStr>, data_len : &mut usize) -> bool {
    let mut i1 : usize;
    let i2 : usize;
    let mut ln_len : usize;
    let wrd_len : usize;
    i1 = file_line.find("#");
    if i1 < MAX_INT {
        *file_line = file_line.substr(0,i1);
    }
    file_line.s = file_line.s.clone() + " ";
    ln_len = file_line.len();
    i1 = file_line.find(":");
    *data_len = 0;
    if i1 < MAX_INT {
        i2 = file_line.find_first_not_of(" -\n\t");
        wrd_len = i1 - i2;
        if headings[0].s == "" || hd_ld_space[0] == i2 {
            headings[0] = file_line.substr(i2,wrd_len);
            hd_ld_space[0] = i2;
            headings[1] = CppStr::from("");
            hd_ld_space[1] = 0;
            headings[2] = CppStr::from("");
            hd_ld_space[2] = 0;
            headings[3] = CppStr::from("");
            hd_ld_space[3] = 0;
        } else if headings[1].s == "" || hd_ld_space[1] == i2 {
            headings[1] = file_line.substr(i2,wrd_len);
            hd_ld_space[1] = i2;
            headings[2] = CppStr::from("");
            hd_ld_space[2] = 0;
            headings[3] = CppStr::from("");
            hd_ld_space[3] = 0;
        } else if headings[2].s == "" || hd_ld_space[2] == i2 {
            headings[2] = file_line.substr(i2,wrd_len);
            hd_ld_space[2] = i2;
            headings[3] = CppStr::from("");
            hd_ld_space[3] = 0;
        } else {
            headings[3] = file_line.substr(i2,wrd_len);
            hd_ld_space[3] = i2;
        }
        i1 += 1usize;
        while i1 < ln_len {
            *file_line = file_line.substr(i1, MAX_INT);
            i1 = file_line.find_first_not_of(" ,[]\t\n");
            if i1 < MAX_INT {
                *file_line = file_line.substr(i1, MAX_INT);
                ln_len = file_line.len();
                i1 = file_line.find_first_of(" ,[]\t\n");
                if i1 < MAX_INT {
                    data[*data_len] = file_line.substr(0,i1);
                    *data_len += 1usize;
                } else {
                    i1 = ln_len;
                }
            } else {
                i1 = ln_len;
            }
        }
        return true;
    } else {
        i1 = file_line.find("- ");
        if i1 < MAX_INT {
            i1 += 1usize;
            while i1 < ln_len {
                *file_line = file_line.substr(i1, MAX_INT);
                i1 = file_line.find_first_not_of(" ,[]\t\n");
                if i1 < MAX_INT {
                    *file_line = file_line.substr(i1, MAX_INT);
                    ln_len = file_line.len();
                    i1 = file_line.find_first_of(" ,[]\t\n");
                    if i1 < MAX_INT {
                        data[*data_len] = file_line.substr(0,i1);
                        *data_len += 1usize;
                    } else {
                        i1 = ln_len;
                    }
                } else {
                    i1 = ln_len;
                }
            }
        }
        return false;
    }
    
}

pub fn print_matrix(mat : &Vec<f64>, size : usize, cols : usize, file_name : &str) {
    let out_file = match File::create(file_name) {
        Err(_why) => panic!("could not open file {}", file_name),
        Ok(file) => file,
    };

    let mut writer = io::BufWriter::new(out_file);

    for i in 0..cols {
        let _ = writer.write(format!("c{},", i).as_bytes());
    }
    let _ = writer.write(b"\n");

    let mut col = 0usize;

    for i in 0..size {
        let _ = writer.write(format!("{},", mat[i]).as_bytes());
        col += 1;
        if col == cols {
            let _ = writer.write(b"\n");
            col = 0;
        }
    }
}