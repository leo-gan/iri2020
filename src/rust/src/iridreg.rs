fn eden_lookup(i: usize, j: usize, k: usize, l: usize, m: usize) -> f32 {
    crate::firi_data::EDEN[i + j * 81 + k * 405 + l * 4860 + m * 53460]
}

pub fn f00(hgt: f32, glat1: f32, iday: i32, zang: f32, f107t: f32) -> Result<(f32, i32), i32> {
    let mut ierror = 0;
    let f107l_check = f107t.max(1.0).min(1000.0).log10();
    
    if hgt < crate::firi_data::TABHE[0] || hgt > crate::firi_data::TABHE[80] ||
       glat1 > crate::firi_data::TABLA[4] || glat1 < -crate::firi_data::TABLA[4] ||
       iday < 1 || iday > 366 ||
       zang < crate::firi_data::TABZA[0] || zang > crate::firi_data::TABZA[10] ||
       f107l_check < crate::firi_data::TABFL[0] || f107l_check > crate::firi_data::TABFL[2] {
        ierror = 2;
    }
    
    // assume height table is in 1 km steps from 60 to 140 km
    let mut i = (hgt as i32) - 59;
    if i < 1 {
        i = 1;
    }
    if i > 80 {
        i = 80;
    }
    let i1 = (i - 1) as usize;
    let i2 = i as usize;
    let h1 = hgt - crate::firi_data::TABHE[i1];
    
    // assume latitude table is in 15 deg steps from 0 to 60 deg
    let mut j = (glat1.abs() as i32) / 15;
    if j < 1 {
        j = 1;
    }
    if j > 4 {
        j = 4;
    }
    let j1 = (j - 1) as usize;
    let j2 = j as usize;
    let deg1 = (glat1.abs() - crate::firi_data::TABLA[j1]) / 15.0;
    
    // assume month table is given for each month
    let tabm = [0, 31, 59, 90, 120, 151, 181, 212, 243, 273, 304, 334];
    let mut mon = 12;
    while mon >= 1 && tabm[mon - 1] > iday {
        mon -= 1;
    }
    let mut day1 = (iday - tabm[mon - 1] - 15) as f32 / 30.0;
    if day1 < 0.0 {
        mon -= 1;
    }
    
    let mut k1;
    let mut k2;
    if mon >= 1 && mon <= 11 {
        k1 = (mon - 1) as usize;
        k2 = mon as usize;
    } else {
        k1 = 11;
        k2 = 0;
    }
    
    // assume zenith angle table has 11 entries between 0 and 130 deg
    let mut l_idx = 1;
    while l_idx < 10 {
        if zang < crate::firi_data::TABZA[l_idx] {
            break;
        }
        l_idx += 1;
    }
    let l1 = l_idx - 1;
    let l2 = l_idx;
    let xhi1 = (zang - crate::firi_data::TABZA[l1]) / (crate::firi_data::TABZA[l2] - crate::firi_data::TABZA[l1]);
    
    // assume solar activity table has 3 entries
    let mut f107l = f107t.max(1.0).min(1000.0).log10();
    f107l = f107l.max(crate::firi_data::TABFL[0]).min(crate::firi_data::TABFL[2]);
    let (m1, m2) = if f107l < crate::firi_data::TABFL[1] {
        (0, 1)
    } else {
        (1, 2)
    };
    let flx1 = (f107l - crate::firi_data::TABFL[m1]) / (crate::firi_data::TABFL[m2] - crate::firi_data::TABFL[m1]);
    
    // Southern hemisphere adjustment
    if glat1 < 0.0 {
        k1 = (k1 + 6) % 12;
        k2 = (k2 + 6) % 12;
    }
    
    let mut edeni = [[[[0.0_f32; 2]; 2]; 2]; 2];
    for (m_i, &m) in [m1, m2].iter().enumerate() {
        for (l_i, &l) in [l1, l2].iter().enumerate() {
            if eden_lookup(i1, j1, k1, l, m) == 0.0 ||
               eden_lookup(i2, j1, k1, l, m) == 0.0 ||
               eden_lookup(i1, j2, k1, l, m) == 0.0 ||
               eden_lookup(i2, j2, k1, l, m) == 0.0 ||
               eden_lookup(i1, j1, k2, l, m) == 0.0 ||
               eden_lookup(i2, j1, k2, l, m) == 0.0 ||
               eden_lookup(i1, j2, k2, l, m) == 0.0 ||
               eden_lookup(i2, j2, k2, l, m) == 0.0 {
                return Err(ierror + 1);
            }
            
            if hgt < crate::firi_data::TABHE[0] {
                edeni[0][0][l_i][m_i] = eden_lookup(i1, j1, k1, l, m);
                edeni[1][0][l_i][m_i] = eden_lookup(i1, j2, k1, l, m);
                edeni[0][1][l_i][m_i] = eden_lookup(i1, j1, k2, l, m);
                edeni[1][1][l_i][m_i] = eden_lookup(i1, j2, k2, l, m);
            } else if hgt > crate::firi_data::TABHE[80] {
                edeni[0][0][l_i][m_i] = eden_lookup(i2, j1, k1, l, m);
                edeni[1][0][l_i][m_i] = eden_lookup(i2, j2, k1, l, m);
                edeni[0][1][l_i][m_i] = eden_lookup(i2, j1, k2, l, m);
                edeni[1][1][l_i][m_i] = eden_lookup(i2, j2, k2, l, m);
            } else {
                edeni[0][0][l_i][m_i] = eden_lookup(i1, j1, k1, l, m) + h1 * (eden_lookup(i2, j1, k1, l, m) - eden_lookup(i1, j1, k1, l, m));
                edeni[1][0][l_i][m_i] = eden_lookup(i1, j2, k1, l, m) + h1 * (eden_lookup(i2, j2, k1, l, m) - eden_lookup(i1, j2, k1, l, m));
                edeni[0][1][l_i][m_i] = eden_lookup(i1, j1, k2, l, m) + h1 * (eden_lookup(i2, j1, k2, l, m) - eden_lookup(i1, j1, k2, l, m));
                edeni[1][1][l_i][m_i] = eden_lookup(i1, j2, k2, l, m) + h1 * (eden_lookup(i2, j2, k2, l, m) - eden_lookup(i1, j2, k2, l, m));
            }
        }
    }
    
    let mut edenij = [[[0.0_f32; 2]; 2]; 2];
    for m_i in 0..2 {
        for l_i in 0..2 {
            if glat1.abs() > crate::firi_data::TABLA[4] {
                edenij[0][l_i][m_i] = edeni[1][0][l_i][m_i];
                edenij[1][l_i][m_i] = edeni[1][1][l_i][m_i];
            } else {
                edenij[0][l_i][m_i] = edeni[0][0][l_i][m_i] + deg1 * (edeni[1][0][l_i][m_i] - edeni[0][0][l_i][m_i]);
                edenij[1][l_i][m_i] = edeni[0][1][l_i][m_i] + deg1 * (edeni[1][1][l_i][m_i] - edeni[0][1][l_i][m_i]);
            }
        }
    }
    
    let mut edenijk = [[0.0_f32; 2]; 2];
    for m_i in 0..2 {
        edenijk[0][m_i] = edenij[0][0][m_i] + day1 * (edenij[1][0][m_i] - edenij[0][0][m_i]);
        edenijk[1][m_i] = edenij[0][1][m_i] + day1 * (edenij[1][1][m_i] - edenij[0][1][m_i]);
    }
    
    let mut edenijkl = [0.0_f32; 2];
    for m_i in 0..2 {
        edenijkl[m_i] = edenijk[0][m_i] + xhi1 * (edenijk[1][m_i] - edenijk[0][m_i]);
    }
    
    let el = edenijkl[0] + flx1 * (edenijkl[1] - edenijkl[0]);
    let edens = 10.0_f32.powf(el);
    
    Ok((edens, ierror))
}
