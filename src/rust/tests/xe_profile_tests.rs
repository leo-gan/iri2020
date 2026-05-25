use iri2020::ioncom::*;
use iri2020::ffi::init_igrf_c;

// Copied aprok from ioncom.rs to print intermediate values
fn debug_aprok(
    j1m: &[i32; 7],
    j2m: &[i32; 7],
    h1: &[[f32; 13]; 7],
    h2: &[[f32; 13]; 7],
    r1m: &[[f32; 13]; 7],
    r2m: &[[f32; 13]; 7],
    rk1m: &[[f32; 13]; 7],
    rk2m: &[[f32; 13]; 7],
    hei: f32,
    xhi: f32,
) -> (f32, f32) {
    let zm = [20.0, 40.0, 60.0, 70.0, 80.0, 85.0, 90.0];
    let h = hei;
    let z = xhi;

    let mut j1 = 1;
    let mut i1 = 0;

    for i in 0..7 {
        i1 = i;
        if (z - zm[i]).abs() < 1e-7 {
            j1 = 0;
        }
        if z <= zm[i] {
            break;
        }
    }

    let eval_at_idx = |idx: usize| -> (f32, f32) {
        let mut i2 = 0;
        let limit1 = j1m[idx] as usize;
        for i in 1..limit1 {
            i2 = i - 1;
            if h < h1[idx][i] {
                break;
            }
            i2 = limit1 - 1;
        }

        let mut i3 = 0;
        let limit2 = j2m[idx] as usize;
        for i in 1..limit2 {
            i3 = i - 1;
            if h < h2[idx][i] {
                break;
            }
            i3 = limit2 - 1;
        }

        let r01 = r1m[idx][i2];
        let r02 = r2m[idx][i3];
        let rk1 = rk1m[idx][i2];
        let rk2 = rk2m[idx][i3];
        let h01 = h1[idx][i2];
        let h02 = h2[idx][i3];

        let r1 = r01 + rk1 * (h - h01);
        let r2 = r02 + rk2 * (h - h02);
        println!("    eval_at_idx: idx={}, i2={}, i3={}, h01={}, h02={}, r01={}, r02={}, rk1={}, rk2={} => r1={}, r2={}", 
                 idx, i2, i3, h01, h02, r01, r02, rk1, rk2, r1, r2);
        (r1, r2)
    };

    if j1 == 0 {
        eval_at_idx(i1)
    } else {
        println!("  aprok interpolating between {} and {}", i1 - 1, i1);
        let (r11, r12) = eval_at_idx(i1);
        let (r1, r2) = eval_at_idx(i1 - 1);
        let rk = (z - zm[i1 - 1]) / (zm[i1] - zm[i1 - 1]);
        let final_r1 = r1 + (r11 - r1) * rk;
        let final_r2 = r2 + (r12 - r2) * rk;
        println!("    rk={}, final_r1={}, final_r2={}", rk, final_r1, final_r2);
        (final_r1, final_r2)
    }
}

fn debug_ionco2(hei: f32, xhi: f32, it: i32, f: f32) -> [f32; 4] {
    let mut z = xhi;
    if z < 20.0 {
        z = 20.0;
    }
    if z > 90.0 {
        z = 90.0;
    }

    let mut r170 = 0.0;
    let mut r270 = 0.0;
    let mut r1140 = 0.0;
    let mut r2140 = 0.0;
    let mut r1 = 0.0;
    let mut r2 = 0.0;

    let is_winter = it == 1 || it == 2 || it == 11 || it == 12;
    let is_summer = it == 5 || it == 6 || it == 7 || it == 8;
    let is_equinox = it == 3 || it == 4 || it == 9 || it == 10;

    if is_winter {
        if f < 140.0 {
            let (v1, v2) = debug_aprok(&J1MW70, &J2MW70, &H1W70, &H2W70, &R1MW70, &R2MW70, &RK1MW70, &RK2MW70, hei, z);
            r170 = v1;
            r270 = v2;
        }
        if f > 70.0 {
            let (v1, v2) = debug_aprok(&J1MW140, &J2MW140, &H1W140, &H2W140, &R1MW140, &R2MW140, &RK1MW140, &RK2MW140, hei, z);
            r1140 = v1;
            r2140 = v2;
        }
        if f <= 70.0 {
            r1 = r170;
            r2 = r270;
        } else if f >= 140.0 {
            r1 = r1140;
            r2 = r2140;
        } else {
            r1 = r170 + (r1140 - r170) * (f - 70.0) / 70.0;
            r2 = r270 + (r2140 - r270) * (f - 70.0) / 70.0;
        }
    } else if is_summer {
        if f < 140.0 {
            let (v1, v2) = debug_aprok(&J1MS70, &J2MS70, &H1S70, &H2S70, &R1MS70, &R2MS70, &RK1MS70, &RK2MS70, hei, z);
            r170 = v1;
            r270 = v2;
        }
        if f > 70.0 {
            let (v1, v2) = debug_aprok(&J1MS140, &J2MS140, &H1S140, &H2S140, &R1MS140, &R2MS140, &RK1MS140, &RK2MS140, hei, z);
            r1140 = v1;
            r2140 = v2;
        }
        if f <= 70.0 {
            r1 = r170;
            r2 = r270;
        } else if f >= 140.0 {
            r1 = r1140;
            r2 = r2140;
        } else {
            r1 = r170 + (r1140 - r170) * (f - 70.0) / 70.0;
            r2 = r270 + (r2140 - r270) * (f - 70.0) / 70.0;
        }
    } else if is_equinox {
        if f < 140.0 {
            println!("Calling debug_aprok for 70:");
            let (v1, v2) = debug_aprok(&J1MR70, &J2MR70, &H1R70, &H2R70, &R1MR70, &R2MR70, &RK1MR70, &RK2MR70, hei, z);
            r170 = v1;
            r270 = v2;
        }
        if f > 70.0 {
            println!("Calling debug_aprok for 140:");
            let (v1, v2) = debug_aprok(&J1MR140, &J2MR140, &H1R140, &H2R140, &R1MR140, &R2MR140, &RK1MR140, &RK2MR140, hei, z);
            r1140 = v1;
            r2140 = v2;
        }
        if f <= 70.0 {
            r1 = r170;
            r2 = r270;
        } else if f >= 140.0 {
            r1 = r1140;
            r2 = r2140;
        } else {
            r1 = r170 + (r1140 - r170) * (f - 70.0) / 70.0;
            r2 = r270 + (r2140 - r270) * (f - 70.0) / 70.0;
        }
    }

    let mut r3 = 0.0;
    let mut r4 = 0.0;
    if hei < 100.0 {
        r3 = 100.0 - (r1 + r2);
    } else {
        r4 = 100.0 - (r1 + r2);
    }
    if r3 < 0.0 {
        r3 = 0.0;
    }
    if r4 < 0.0 {
        r4 = 0.0;
    }

    println!("    Before rounding: r1={}, r2={}, r3={}, r4={}", r1, r2, r3, r4);

    [
        r1.round(),
        r2.round(),
        r3.round(),
        r4.round(),
    ]
}

#[test]
fn test_ioncom_equivalence() {
    unsafe {
        init_igrf_c();
    }

    let hx = 150.0;
    let zd = 75.0;
    let fs = 70.0;
    let ismo = 4;
    let id = 15;
    let fd = 0.0;

    let molecular = debug_ionco2(hx, zd, ismo, fs);
    println!("DEBUG: final molecular output = {:?}", molecular);

    let mut dion = [0.0_f32; 7];
    unsafe {
        use iri2020::ffi::iondani_c;
        iondani_c(id, ismo, hx, zd, fd, fs, dion.as_mut_ptr());
    }
    println!("DEBUG: Fortran iondani output = {:?}", dion);
    println!("DEBUG: Fortran rounded: O+={}, NO+={}, O2+={}, Cluster+={}", 
        dion[0].round(), dion[4].round(), dion[5].round(), dion[6].round());

    assert_eq!(molecular[0], dion[4].round(), "NO+ mismatch");
    assert_eq!(molecular[1], dion[5].round(), "O2+ mismatch");
    assert_eq!(molecular[2], dion[6].round(), "Cluster+ mismatch");
    assert_eq!(molecular[3], dion[0].round(), "O+ mismatch");
}

#[test]
fn test_xe_profile_equivalence() {
    use iri2020::xe_profile::XeProfile;
    use iri2020::ffi::{set_xe_blocks_c, xe_1_c};

    unsafe {
        init_igrf_c();
    }

    let profile = XeProfile {
        hmf2: 300.0,
        xnmf2: 1.0e12,
        hmf1: 200.0,
        f1reg: true,
        b0: 80.0,
        b1: 2.0,
        c1: 0.1,
        hz: 180.0,
        t: 0.0,
        hst: 120.0,
        hme: 110.0,
        xnme: 1.0e11,
        hef: 115.0,
        night: false,
        e: [0.01, 0.02, 0.03, 0.04],
        hmd: 80.0,
        xnmd: 1.0e9,
        hdx: 95.0,
        d1: 0.05,
        xkk: 0.6,
        fp30: 0.1,
        fp3u: 0.2,
        fp1: 0.3,
        fp2: 0.4,
        beta: 1.2,
        eta: 0.8,
        delta: 2.5,
        zeta: 1.5,
        b2top: 140.0,
        itopn: 0,
        tcor1: 0.0,
        tcor2: 0.0,
    };

    unsafe {
        set_xe_blocks_c(
            profile.hmf2,
            profile.xnmf2,
            profile.hmf1,
            profile.f1reg,
            profile.b0,
            profile.b1,
            profile.c1,
            profile.hz,
            profile.t,
            profile.hst,
            profile.hme,
            profile.xnme,
            profile.hef,
            profile.night,
            profile.e.as_ptr(),
            profile.hmd,
            profile.xnmd,
            profile.hdx,
            profile.d1,
            profile.xkk,
            profile.fp30,
            profile.fp3u,
            profile.fp1,
            profile.fp2,
            profile.beta,
            profile.eta,
            profile.delta,
            profile.zeta,
            profile.b2top,
            profile.itopn,
            profile.tcor1,
            profile.tcor2,
        );
    }

    let mut hmf2_f = 0.0;
    let mut xnmf2_f = 0.0;
    let mut hmf1_f = 0.0;
    let mut f1reg_f = false;
    let mut b0_f = 0.0;
    let mut b1_f = 0.0;
    let mut c1_f = 0.0;
    let mut hz_f = 0.0;
    let mut t_f = 0.0;
    let mut hst_f = 0.0;
    let mut hme_f = 0.0;
    let mut xnme_f = 0.0;
    let mut hef_f = 0.0;
    let mut night_f = false;
    let mut e_f = [0.0; 4];
    let mut hmd_f = 0.0;
    let mut xnmd_f = 0.0;
    let mut hdx_f = 0.0;
    let mut d1_f = 0.0;
    let mut xkk_f = 0.0;
    let mut fp30_f = 0.0;
    let mut fp3u_f = 0.0;
    let mut fp1_f = 0.0;
    let mut fp2_f = 0.0;
    let mut beta_f = 0.0;
    let mut eta_f = 0.0;
    let mut delta_f = 0.0;
    let mut zeta_f = 0.0;
    let mut b2top_f = 0.0;
    let mut itopn_f = 0;
    let mut tcor1_f = 0.0;
    let mut tcor2_f = 0.0;

    unsafe {
        use iri2020::ffi::get_xe_blocks_c;
        get_xe_blocks_c(
            &mut hmf2_f,
            &mut xnmf2_f,
            &mut hmf1_f,
            &mut f1reg_f,
            &mut b0_f,
            &mut b1_f,
            &mut c1_f,
            &mut hz_f,
            &mut t_f,
            &mut hst_f,
            &mut hme_f,
            &mut xnme_f,
            &mut hef_f,
            &mut night_f,
            e_f.as_mut_ptr(),
            &mut hmd_f,
            &mut xnmd_f,
            &mut hdx_f,
            &mut d1_f,
            &mut xkk_f,
            &mut fp30_f,
            &mut fp3u_f,
            &mut fp1_f,
            &mut fp2_f,
            &mut beta_f,
            &mut eta_f,
            &mut delta_f,
            &mut zeta_f,
            &mut b2top_f,
            &mut itopn_f,
            &mut tcor1_f,
            &mut tcor2_f,
        );
    }
    println!("Fortran readback: hmf2={}, xnmf2={}, hmf1={}, f1reg={}, b0={}, b1={}, c1={}, hz={}, t={}, hst={}, hme={}, xnme={}, hef={}, night={}",
        hmf2_f, xnmf2_f, hmf1_f, f1reg_f, b0_f, b1_f, c1_f, hz_f, t_f, hst_f, hme_f, xnme_f, hef_f, night_f
    );

    let h = 117.0;
    unsafe {
        use iri2020::ffi::{xe1_c, xe2_c, xe3_1_c, xe4_1_c, xe5_c, xe6_c};
        println!("For h = 117.0:");
        println!("  xe_1:   rust = {:e}, fortran = {:e}", profile.xe_1(h), xe_1_c(h));
        println!("  xe1:    rust = {:e}, fortran = {:e}", profile.xe1(h), xe1_c(h));
        println!("  xe2:    rust = {:e}, fortran = {:e}", profile.xe2(h), xe2_c(h));
        println!("  xe3_1:  rust = {:e}, fortran = {:e}", profile.xe3_1(h), xe3_1_c(h));
        println!("  xe4_1:  rust = {:e}, fortran = {:e}", profile.xe4_1(h), xe4_1_c(h));
        println!("  xe5:    rust = {:e}, fortran = {:e}", profile.xe5(h), xe5_c(h));
        println!("  xe6:    rust = {:e}, fortran = {:e}", profile.xe6(h), xe6_c(h));
    }

    // Let's test a range of heights
    let heights = [
        80.0, 90.0, 100.0, 110.0, 112.0, 117.0, 120.0, 130.0, 150.0, 190.0, 200.0, 220.0, 250.0,
        300.0, 350.0, 500.0, 1000.0,
    ];
    for &h in &heights {
        let rust_val = profile.xe_1(h);
        let fortran_val = unsafe { xe_1_c(h) };
        println!("h = {}, rust = {:e}, fortran = {:e}", h, rust_val, fortran_val);
        if fortran_val.is_infinite() && rust_val.is_infinite() {
            // Both are infinite, which is expected for these dummy inputs
        } else if fortran_val != 0.0 {
            let diff = (rust_val - fortran_val).abs() / fortran_val;
            assert!(
                diff < 1e-4,
                "h = {}, rust = {:e}, fortran = {:e}, relative diff = {:e}",
                h,
                rust_val,
                fortran_val,
                diff
            );
        } else {
            assert_eq!(rust_val, 0.0);
        }
    }
}


