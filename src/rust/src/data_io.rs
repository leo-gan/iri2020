use std::fs::File;
use std::io::{self, BufRead, BufReader};
use std::path::Path;

pub struct IgRzData {
    pub aig: [f32; 806],
    pub arz: [f32; 806],
    pub iymst: i32,
    pub iymend: i32,
}

impl IgRzData {
    pub fn load(data_dir: &str) -> io::Result<Self> {
        let path = Path::new(data_dir).join("ig_rz.dat");
        let content = std::fs::read_to_string(path)?;
        
        let tokens: Vec<&str> = content
            .split(|c: char| c == ',' || c == '\n' || c == '\r' || c.is_whitespace())
            .map(|s| s.trim())
            .filter(|s| !s.is_empty())
            .collect();
            
        if tokens.len() < 7 {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                "ig_rz.dat does not contain enough header elements",
            ));
        }
        
        let _iupd = tokens[0].parse::<i32>().map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;
        let iupm = tokens[1].parse::<i32>().map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;
        let iupy = tokens[2].parse::<i32>().map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;
        
        let imst = tokens[3].parse::<i32>().map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;
        let iyst = tokens[4].parse::<i32>().map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;
        let imend = tokens[5].parse::<i32>().map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;
        let iyend = tokens[6].parse::<i32>().map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;
        
        let iymst = iyst * 100 + imst;
        let iymend = iyend * 100 + imend;
        
        let inum_vals = 3 - imst + (iyend - iyst) * 12 + imend;
        let inum_vals_usize = inum_vals as usize;
        
        if tokens.len() < 7 + 2 * inum_vals_usize {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                format!(
                    "ig_rz.dat does not contain enough data elements: expected {}, got {}",
                    7 + 2 * inum_vals_usize,
                    tokens.len()
                ),
            ));
        }
        
        let mut aig = [0.0f32; 806];
        let mut arz = [0.0f32; 806];
        
        for i in 0..inum_vals_usize {
            aig[i] = tokens[7 + i].parse::<f32>().map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;
            arz[i] = tokens[7 + inum_vals_usize + i].parse::<f32>().map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;
        }
        
        if iupy * 100 + iupm > 201609 {
            let inum_chan = 3 - imst + (2014 - iyst) * 12;
            let inum_chan_usize = (inum_chan as usize - 1).max(0);
            for jj in inum_chan_usize..inum_vals_usize {
                arz[jj] *= 0.7;
            }
        }
        
        Ok(Self { aig, arz, iymst, iymend })
    }
}

pub struct Apf107Data {
    pub aap: Vec<[i32; 9]>,
    pub af107: Vec<[f32; 3]>,
    pub n: i32,
}

impl Apf107Data {
    pub fn load(data_dir: &str) -> io::Result<Self> {
        let path = Path::new(data_dir).join("apf107.dat");
        let file = File::open(path)?;
        let reader = BufReader::new(file);
        
        let mut aap = vec![[0i32; 9]; 27000];
        let mut af107 = vec![[0.0f32; 3]; 27000];
        let mut n = 0;
        
        for line in reader.lines() {
            let line = line?;
            if line.trim().is_empty() {
                continue;
            }
            if let Some((aap_row, af_row)) = parse_apf107_line(&line) {
                if n < 27000 {
                    aap[n] = aap_row;
                    af107[n] = af_row;
                    n += 1;
                }
            }
        }
        
        Ok(Self { aap, af107, n: n as i32 })
    }
}

fn parse_apf107_line(line: &str) -> Option<([i32; 9], [f32; 3])> {
    let mut padded = line.to_string();
    if padded.len() < 54 {
        padded.push_str(&" ".repeat(54 - padded.len()));
    }
    
    let mut aap = [0i32; 9];
    for j in 0..8 {
        let start = 9 + j * 3;
        let s = padded[start..start+3].trim();
        aap[j] = s.parse::<i32>().unwrap_or(0);
    }
    let s_iapda = padded[33..36].trim();
    aap[8] = s_iapda.parse::<i32>().unwrap_or(0);
    
    let f107d = padded[39..44].trim().parse::<f32>().unwrap_or(0.0);
    let mut f107_81 = padded[44..49].trim().parse::<f32>().unwrap_or(0.0);
    let mut f107_365 = padded[49..54].trim().parse::<f32>().unwrap_or(0.0);
    
    if f107_81 < -4.0 {
        f107_81 = f107d;
    }
    if f107_365 < -4.0 {
        f107_365 = f107d;
    }
    
    Some((aap, [f107d, f107_81, f107_365]))
}

pub struct McsatData;

impl McsatData {
    pub fn load(data_dir: &str, month: i32) -> io::Result<[[f64; 48]; 149]> {
        let month_str = format!("{:02}", month + 10);
        let filename = format!("mcsat{}.dat", month_str);
        let path = Path::new(data_dir).join(filename);
        let content = std::fs::read_to_string(path)?;
        
        let mut coeff = [[0.0f64; 48]; 149];
        
        let lines: Vec<&str> = content.lines().collect();
        if lines.len() < 1200 {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                format!("mcsat file does not have enough lines: expected 1200, got {}", lines.len()),
            ));
        }
        
        for j in 0..48 {
            let mut values = Vec::with_capacity(149);
            for line_idx in 0..25 {
                let actual_line_idx = j * 25 + line_idx;
                let line = lines[actual_line_idx];
                
                for k in 0..6 {
                    let start = k * 12;
                    if start >= line.len() {
                        break;
                    }
                    let end = (start + 12).min(line.len());
                    let val_str = line[start..end].trim();
                    if val_str.is_empty() {
                        continue;
                    }
                    let val_clean = val_str.replace('D', "E").replace('d', "e");
                    if let Ok(val) = val_clean.parse::<f64>() {
                        values.push(val);
                    }
                }
            }
            if values.len() < 149 {
                return Err(io::Error::new(
                    io::ErrorKind::InvalidData,
                    format!("Hour {} has only {} values, expected 149", j, values.len()),
                ));
            }
            for i in 0..149 {
                coeff[i][j] = values[i];
            }
        }
        
        Ok(coeff)
    }
}

pub struct CcirUrsiData;

impl CcirUrsiData {
    pub fn load(
        data_dir: &str,
        month: i32,
        is_ccir: bool,
    ) -> io::Result<(
        [[[f32; 2]; 76]; 13],
        Option<[[[f32; 2]; 49]; 9]>,
    )> {
        let month_str = format!("{:02}", month + 10);
        let prefix = if is_ccir { "ccir" } else { "ursi" };
        let filename = format!("{}{}.asc", prefix, month_str);
        let path = Path::new(data_dir).join(filename);
        let content = std::fs::read_to_string(path)?;
        
        let mut flat_vals = Vec::new();
        for line in content.lines() {
            flat_vals.extend(parse_asc_line(line));
        }
        
        if flat_vals.len() < 1976 {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                format!("Not enough floats in file: expected at least 1976, got {}", flat_vals.len()),
            ));
        }
        
        let mut f2 = [[[0.0f32; 2]; 76]; 13];
        let mut idx = 0;
        for k in 0..2 {
            for j in 0..76 {
                for i in 0..13 {
                    f2[i][j][k] = flat_vals[idx];
                    idx += 1;
                }
            }
        }
        
        let mut fm3 = None;
        if is_ccir {
            if flat_vals.len() < 1976 + 882 {
                return Err(io::Error::new(
                    io::ErrorKind::InvalidData,
                    format!("Not enough floats in CCIR file: expected 2858, got {}", flat_vals.len()),
                ));
            }
            let mut fm3_arr = [[[0.0f32; 2]; 49]; 9];
            for k in 0..2 {
                for j in 0..49 {
                    for i in 0..9 {
                        fm3_arr[i][j][k] = flat_vals[idx];
                        idx += 1;
                    }
                }
            }
            fm3 = Some(fm3_arr);
        }
        
        Ok((f2, fm3))
    }
}

fn parse_asc_line(line: &str) -> Vec<f32> {
    let mut values = Vec::new();
    if line.len() <= 1 {
        return values;
    }
    let s = &line[1..];
    let chunk_size = 15;
    for i in 0..4 {
        let start = i * chunk_size;
        if start >= s.len() {
            break;
        }
        let end = (start + chunk_size).min(s.len());
        let chunk = &s[start..end];
        let chunk_trimmed = chunk.trim();
        if chunk_trimmed.is_empty() {
            continue;
        }
        let chunk_clean = chunk_trimmed.replace('D', "E").replace('d', "e");
        if let Ok(val) = chunk_clean.parse::<f32>() {
            values.push(val);
        }
    }
    values
}
