use std::os::raw::{c_int, c_float};

extern "C" {
    pub fn read_ig_rz_c();
    pub fn readapf107_c();

    pub fn iri_sub_c(
        jf_c: *const bool,
        jmag_c: c_int,
        alati_c: c_float,
        along_c: c_float,
        iyyyy_c: c_int,
        mmdd_c: c_int,
        dhour_c: c_float,
        heibeg_c: c_float,
        heiend_c: c_float,
        heistp_c: c_float,
        outf_c: *mut c_float,
        oarr_c: *mut c_float,
    );

    pub fn iritec_c(
        alati_c: c_float,
        along_c: c_float,
        jmag_c: c_int,
        jf_c: *const bool,
        iyyyy_c: c_int,
        mmdd_c: c_int,
        hour_c: c_float,
        hbgn_c: c_float,
        hend_c: c_float,
        hstep_c: c_float,
        oarr_c: *mut c_float,
        tec_c: *mut c_float,
        tect_c: *mut c_float,
    );

    pub fn feldcof_c(year_c: c_float);

    pub fn feldg_c(
        glat_c: c_float,
        glon_c: c_float,
        alt_c: c_float,
        bnorth_c: *mut c_float,
        beast_c: *mut c_float,
        bdown_c: *mut c_float,
        babs_c: *mut c_float,
    );

    pub fn igrf_c(
        iy_c: c_int,
        nm_c: c_int,
        r_c: c_float,
        t_c: c_float,
        f_c: c_float,
        br_c: *mut c_float,
        bt_c: *mut c_float,
        bf_c: *mut c_float,
    );

    pub fn igrf_dip_c(
        xlat_c: c_float,
        xlong_c: c_float,
        year_c: c_float,
        height_c: c_float,
        dec_c: *mut c_float,
        dip_c: *mut c_float,
        dipl_c: *mut c_float,
        ymodip_c: *mut c_float,
    );

    pub fn init_igrf_c();

    pub fn gtd7_c(
        iyd_c: c_int,
        sec_c: c_float,
        alt_c: c_float,
        glat_c: c_float,
        glong_c: c_float,
        stl_c: c_float,
        f107a_c: c_float,
        f107_c: c_float,
        ap_c: *const c_float,
        mass_c: c_int,
        d_c: *mut c_float,
        t_c: *mut c_float,
    );

    pub fn gtd7d_c(
        iyd_c: c_int,
        sec_c: c_float,
        alt_c: c_float,
        glat_c: c_float,
        glong_c: c_float,
        stl_c: c_float,
        f107a_c: c_float,
        f107_c: c_float,
        ap_c: *const c_float,
        mass_c: c_int,
        d_c: *mut c_float,
        t_c: *mut c_float,
    );

    pub fn tselec_c(sv_c: *const c_float);

    pub fn meters_c(meter_c: bool);

    pub fn vfjmodelrocstart_c(vzm_c: *mut c_float);

    pub fn vfjmodelrocinit_c(
        f107_c: c_float,
        idoy_c: c_int,
        jseas_c: *mut c_int,
        jsfl_c: *mut c_int,
    );

    pub fn vfjmodelroc_c(
        fjm_c: *const c_float,
        ttl_c: c_float,
        gglon_c: c_float,
        jseas_c: c_int,
        jsfl_c: c_int,
        viv_c: *mut c_float,
    );

    pub fn get_igrz_c(
        aig_c: *mut c_float,
        arz_c: *mut c_float,
        iymst_c: *mut c_int,
        iymend_c: *mut c_int,
    );

    pub fn get_apfa_c(
        aap_c: *mut c_int,
        af107_c: *mut c_float,
        n_c: *mut c_int,
    );

    pub fn read_data_sd_c(
        month_c: c_int,
        coeff_month_c: *mut f64,
    );

    pub fn read_coeff_c(
        month_c: c_int,
        is_ccir_c: bool,
        f2_c: *mut c_float,
        fm3_c: *mut c_float,
    );
}
