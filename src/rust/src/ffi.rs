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

    pub fn shellg_c(
        glat_c: c_float,
        glon_c: c_float,
        alt_c: c_float,
        fl_c: *mut c_float,
        icode_c: *mut c_int,
        b0_c: *mut c_float,
    );

    pub fn iondani_c(
        id_c: c_int,
        ismo_c: c_int,
        hx_c: c_float,
        zd_c: c_float,
        fd_c: c_float,
        fs_c: c_float,
        dion_c: *mut c_float,
    );

    pub fn f1_prob_c(
        sza_c: c_float,
        glat_c: c_float,
        rz12_c: c_float,
        f1prob_c: *mut c_float,
        f1probl_c: *mut c_float,
    );

    pub fn set_xe_blocks_c(
        hmf2_c: c_float,
        xnmf2_c: c_float,
        hmf1_c: c_float,
        f1reg_c: bool,
        b0_c: c_float,
        b1_c: c_float,
        c1_c: c_float,
        hz_c: c_float,
        t_c: c_float,
        hst_c: c_float,
        hme_c: c_float,
        xnme_c: c_float,
        hef_c: c_float,
        night_c: bool,
        e_c: *const c_float,
        hmd_c: c_float,
        xnmd_c: c_float,
        hdx_c: c_float,
        d1_c: c_float,
        xkk_c: c_float,
        fp30_c: c_float,
        fp3u_c: c_float,
        fp1_c: c_float,
        fp2_c: c_float,
        beta_c: c_float,
        eta_c: c_float,
        delta_c: c_float,
        zeta_c: c_float,
        b2top_c: c_float,
        itopn_c: c_int,
        tcor1_c: c_float,
        tcor2_c: c_float,
    );

    pub fn get_xe_blocks_c(
        hmf2_c: *mut c_float,
        xnmf2_c: *mut c_float,
        hmf1_c: *mut c_float,
        f1reg_c: *mut bool,
        b0_c: *mut c_float,
        b1_c: *mut c_float,
        c1_c: *mut c_float,
        hz_c: *mut c_float,
        t_c: *mut c_float,
        hst_c: *mut c_float,
        hme_c: *mut c_float,
        xnme_c: *mut c_float,
        hef_c: *mut c_float,
        night_c: *mut bool,
        e_c: *mut c_float,
        hmd_c: *mut c_float,
        xnmd_c: *mut c_float,
        hdx_c: *mut c_float,
        d1_c: *mut c_float,
        xkk_c: *mut c_float,
        fp30_c: *mut c_float,
        fp3u_c: *mut c_float,
        fp1_c: *mut c_float,
        fp2_c: *mut c_float,
        beta_c: *mut c_float,
        eta_c: *mut c_float,
        delta_c: *mut c_float,
        zeta_c: *mut c_float,
        b2top_c: *mut c_float,
        itopn_c: *mut c_int,
        tcor1_c: *mut c_float,
        tcor2_c: *mut c_float,
    );

    pub fn xe_1_c(h_c: c_float) -> c_float;
    pub fn xe1_c(h_c: c_float) -> c_float;
    pub fn xe2_c(h_c: c_float) -> c_float;
    pub fn xe3_1_c(h_c: c_float) -> c_float;
    pub fn xe4_1_c(h_c: c_float) -> c_float;
    pub fn xe5_c(h_c: c_float) -> c_float;
    pub fn xe6_c(h_c: c_float) -> c_float;
    pub fn dxe1n_c(h_c: c_float) -> c_float;
    pub fn topq_c(h_c: c_float, no_c: c_float, hmax_c: c_float, ho_c: c_float) -> c_float;
    pub fn zero_c(delta_c: c_float) -> c_float;

    pub fn foeedi_c(cov_c: c_float, xhi_c: c_float, xhim_c: c_float, xlati_c: c_float) -> c_float;
    pub fn xmded_c(xhi_c: c_float, r_c: c_float, yw_c: c_float) -> c_float;
    pub fn valgul_c(
        xhi_c: c_float,
        hvb_c: *mut c_float,
        vwu_c: *mut c_float,
        vwa_c: *mut c_float,
        vdp_c: *mut c_float,
    );

    pub fn spharm_c(c_c: *mut c_float, l_c: c_int, m_c: c_int, colat_c: c_float, az_c: c_float);
    pub fn spharm_ik_c(c_c: *mut c_float, l_c: c_int, m_c: c_int, colat_c: c_float, az_c: c_float);

    pub fn dregion_c(
        z_c: c_float,
        it_c: c_int,
        f_c: c_float,
        vkp_c: c_float,
        f5sw_c: c_float,
        f6wa_c: c_float,
        elg_c: *mut c_float,
    );

    pub fn f00_c(
        hgt_c: c_float,
        glat1_c: c_float,
        iday_c: c_int,
        zang_c: c_float,
        f107t_c: c_float,
        edens_c: *mut c_float,
        ierror_c: *mut c_int,
    );

    pub fn shamdb0d_c(rlat_c: c_float, flon_c: c_float, t_c: c_float, rz_c: c_float) -> c_float;
    pub fn shab1d_c(flat_c: c_float, flon_c: c_float, t_c: c_float, rz_c: c_float) -> c_float;

    pub fn geocgm01_c(
        icor_c: c_int,
        iyear_c: c_int,
        hi_c: c_float,
        dat_c: *mut c_float,
        pla_c: *mut c_float,
        plo_c: *mut c_float,
    );
}

