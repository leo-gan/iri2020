pub fn spharm_ik(c: &mut [f32], l: i32, m: i32, colat: f32, az: f32) {
    c[0] = 1.0;
    let mut k = 2;
    let x = colat.cos();
    c[k - 1] = x;
    k += 1;
    for i in 2..=l {
        c[k - 1] = ((2 * i - 1) as f32 * x * c[k - 2] - (i - 1) as f32 * c[k - 3]) / i as f32;
        k += 1;
    }
    let y = colat.sin();
    for mt in 1..=m {
        let caz = (mt as f32 * az).cos();
        let saz = (mt as f32 * az).sin();
        c[k - 1] = y.powi(mt);
        k += 1;
        if mt == l {
            let n = (l - mt + 1) as usize;
            for _ in 1..=n {
                let temp = c[k - 1 - n];
                c[k - 1] = temp * saz;
                c[k - 1 - n] = temp * caz;
                k += 1;
            }
            continue;
        }
        c[k - 1] = c[k - 2] * x * (2 * mt + 1) as f32;
        k += 1;
        if (mt + 1) == l {
            let n = (l - mt + 1) as usize;
            for _ in 1..=n {
                let temp = c[k - 1 - n];
                c[k - 1] = temp * saz;
                c[k - 1 - n] = temp * caz;
                k += 1;
            }
            continue;
        }
        for i in (2 + mt)..=l {
            c[k - 1] = ((2 * i - 1) as f32 * x * c[k - 2] - (i + mt - 1) as f32 * c[k - 3]) / (i - mt) as f32;
            k += 1;
        }
        let n = (l - mt + 1) as usize;
        for _ in 1..=n {
            let temp = c[k - 1 - n];
            c[k - 1] = temp * saz;
            c[k - 1 - n] = temp * caz;
            k += 1;
        }
    }
}

pub fn spharm(c: &mut [f32], l: i32, m: i32, colat: f32, az: f32) {
    c[0] = 1.0;
    let mut k = 2;
    let x = colat.cos();
    c[k - 1] = x;
    k += 1;
    for i in 2..=l {
        c[k - 1] = ((2 * i - 1) as f32 * x * c[k - 2] - (i - 1) as f32 * c[k - 3]) / i as f32;
        k += 1;
    }
    let y = colat.sin();
    for mt in 1..=m {
        let caz = (mt as f32 * az).cos();
        let saz = (mt as f32 * az).sin();
        c[k - 1] = y.powi(mt);
        k += 1;
        if mt == l {
            let n = (l - mt + 1) as usize;
            for _ in 1..=n {
                let temp = c[k - 1 - n];
                c[k - 1] = temp * caz;
                c[k - 1 - n] = temp * saz;
                k += 1;
            }
            continue;
        }
        c[k - 1] = c[k - 2] * x * (2 * mt + 1) as f32;
        k += 1;
        if (mt + 1) == l {
            let n = (l - mt + 1) as usize;
            for _ in 1..=n {
                let temp = c[k - 1 - n];
                c[k - 1] = temp * caz;
                c[k - 1 - n] = temp * saz;
                k += 1;
            }
            continue;
        }
        for i in (2 + mt)..=l {
            c[k - 1] = ((2 * i - 1) as f32 * x * c[k - 2] - (i + mt - 1) as f32 * c[k - 3]) / (i - mt) as f32;
            k += 1;
        }
        let n = (l - mt + 1) as usize;
        for _ in 1..=n {
            let temp = c[k - 1 - n];
            c[k - 1] = temp * caz;
            c[k - 1 - n] = temp * saz;
            k += 1;
        }
    }
}
