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
}
