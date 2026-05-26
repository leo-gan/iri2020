pub fn ioncorr(tec: f32, f: f32) -> f32 {
    if f == 0.0 {
        0.0
    } else {
        40.3 * tec / (f * f)
    }
}

pub fn iritec(
    alati: f32,
    along: f32,
    jmag: i32,
    jf: &[bool; 50],
    iy: i32,
    md: i32,
    hour: f32,
    hbeg: f32,
    hend: f32,
    hstep: f32,
    oarr: &mut [f32; 100],
    tecbo: &mut f32,
    tecto: &mut f32,
) {
    let mut jff = *jf;
    jff[1] = false;  // jf(2) = .false. (index 1)
    jff[2] = false;  // jf(3) = .false. (index 2)
    jff[20] = false; // jf(21) = .false. (index 20)
    jff[27] = false; // jf(28) = .false. (index 27)
    jff[32] = false; // jf(33) = .false. (index 32)
    jff[33] = false; // jf(34) = .false. (index 33)
    jff[34] = false; // jf(35) = .false. (index 34)
    jff[46] = false; // jf(47) = .false. (index 46)

    let iisect = (((hend - hbeg) / hstep) / 1000.0) as i32;
    let hsect = 1000.0 * hstep;
    let hastep = hstep / 2.0;
    let mut tect = 0.0_f32;
    let mut tecb = 0.0_f32;

    let mut hm_f2 = 0.0_f32;
    let mut xnm_f2 = 0.0_f32;
    let mut xnorm = 1.0_f32;

    let mut outf = [0.0_f32; 20 * 1000];

    if iisect >= 1 {
        for j in 1..=iisect {
            let abeg = hbeg + (j - 1) as f32 * hsect + hastep;
            let aend = abeg + hsect - hstep;

            unsafe {
                crate::ffi::iri_sub_c(
                    jff.as_ptr(),
                    jmag,
                    alati,
                    along,
                    iy,
                    md,
                    hour,
                    abeg,
                    aend,
                    hstep,
                    outf.as_mut_ptr(),
                    oarr.as_mut_ptr(),
                );
            }

            if j == 1 {
                hm_f2 = oarr[1];
                xnm_f2 = oarr[0];
                xnorm = if xnm_f2 == 0.0 { 1e-10 } else { xnm_f2 / 1000.0 };
            }

            let mut hx = abeg + hastep;
            for jj in 1..=1000 {
                let yyy = outf[((jj - 1) * 20) as usize] * hstep / xnorm;
                if hx <= hm_f2 {
                    tecb += yyy;
                } else {
                    tect += yyy;
                }
                hx += hstep;
            }
        }
    }

    // Last segment
    let hlastbeg = hbeg + iisect as f32 * hsect;
    let ilast = ((hend - hlastbeg) / hstep) as i32;

    let abeg = hlastbeg + hastep;
    let aend = hlastbeg + ilast as f32 * hstep - hastep;

    unsafe {
        crate::ffi::iri_sub_c(
            jf.as_ptr(),
            jmag,
            alati,
            along,
            iy,
            md,
            hour,
            abeg,
            aend,
            hstep,
            outf.as_mut_ptr(),
            oarr.as_mut_ptr(),
        );
    }

    if iisect < 1 {
        hm_f2 = oarr[1];
        xnm_f2 = oarr[0];
        xnorm = if xnm_f2 == 0.0 { 1e-10 } else { xnm_f2 / 1000.0 };
    }

    let mut hx = abeg + hastep;
    for jj in 1..=ilast {
        let yyy = outf[((jj - 1) * 20) as usize] * hstep / xnorm;
        if hx <= hm_f2 {
            tecb += yyy;
        } else {
            tect += yyy;
        }
        hx += hstep;
    }

    *tecto = tect * xnm_f2;
    *tecbo = tecb * xnm_f2;
}
