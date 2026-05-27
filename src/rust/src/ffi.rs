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

    pub fn fout_c(
        modip_c: c_float,
        lati_c: c_float,
        longi_c: c_float,
        hourut_c: c_float,
        ff0_c: *const c_float,
    ) -> c_float;

    pub fn xmout_c(
        modip_c: c_float,
        lati_c: c_float,
        longi_c: c_float,
        hourut_c: c_float,
        xm0_c: *const c_float,
    ) -> c_float;

    pub fn apf_c(
        isdate_c: c_int,
        hourut_c: c_float,
        indap_c: *mut c_int,
    );

    pub fn apfmsis_c(
        isdate_c: c_int,
        hourut_c: c_float,
        iapo_c: *mut c_float,
    );

    pub fn apf_only_c(
        iyyyy_c: c_int,
        mm_c: c_int,
        id_c: c_int,
        f107d_c: *mut c_float,
        f107pd_c: *mut c_float,
        f10781_c: *mut c_float,
        f107365_c: *mut c_float,
        iapd_c: *mut c_int,
        isdate_c: *mut c_int,
    );

    pub fn storm_c(
        indap_c: *const c_int,
        lati_c: c_float,
        longi_c: c_float,
        icoord_c: c_int,
        cglat_c: *mut c_float,
        kut_c: c_int,
        daynr_c: c_int,
        stormcorr_c: *mut c_float,
    );

    pub fn storme_ap_c(
        daynr_c: c_int,
        mlat_c: c_float,
        ah3_c: c_float,
    ) -> c_float;

    pub fn auroral_boundary_c(
        xkp_c: c_float,
        xmlt_c: c_float,
        cgmlat_c: *mut c_float,
        ab_mlat_c: *mut c_float,
    );

    pub fn shamdhmf2_c(
        rlat_c: c_float,
        flon_c: c_float,
        t_c: c_float,
        rz_c: c_float,
        hmf2_c: *mut c_float,
    );

    pub fn model_hmf2_c(
        iday_c: c_int,
        month_c: c_int,
        ut_c: c_float,
        modip_c: c_float,
        longi_c: c_float,
        f10781_c: c_float,
        hmf2_c: *mut c_float,
    );

    pub fn hmf2ed_c(
        magbr_c: c_float,
        rssn_c: c_float,
        ratf_c: c_float,
        xm3000_c: c_float,
    ) -> c_float;

    pub fn toph05_c(
        cov_c: c_float,
        mlat_c: c_float,
        hour_c: c_float,
        hmf2_c: c_float,
        hei05_c: *mut c_float,
        sday_c: c_float,
    );

    pub fn elteik_c(
        pf107y_c: c_int,
        invdip_c: c_float,
        mlt_c: c_float,
        ddd_c: c_float,
        pf107_c: c_float,
        tev_c: *mut c_float,
        sigtev_c: *mut c_float,
    );

    pub fn iontif_c(
        pf107y_c: c_int,
        invdip_c: c_float,
        mlt_c: c_float,
        ddd_c: c_float,
        pf107_c: c_float,
        tiv_c: *mut c_float,
        sigtiv_c: *mut c_float,
    );

    pub fn gallden_c(
        l_c: c_float,
        day_c: c_float,
        rz12_c: c_float,
    ) -> c_float;

    pub fn ohzden_c(
        l_c: c_float,
        lat_c: c_float,
    ) -> c_float;

    pub fn fof1ed_c(
        absmlt_c: c_float,
        rssn_c: c_float,
        xhi3_c: c_float,
    ) -> c_float;

    pub fn f1_c1_c(
        absmdp_c: c_float,
        hour_c: c_float,
        sax2_c: c_float,
        sux2_c: c_float,
    ) -> c_float;

    pub fn tal_c(
        hdeep_c: c_float,
        depth_c: c_float,
        width_c: c_float,
        dlndh_c: c_float,
        ext_c: *mut bool,
        e_c: *mut c_float,
    );

    pub fn rogul_c(
        day_c: c_float,
        xhi_c: c_float,
        seax_c: *mut c_float,
        grat_c: *mut c_float,
    );

    pub fn inilay_c(
        night_c: bool,
        f1reg_c: bool,
        nmf2_c: c_float,
        nmf1_c: c_float,
        nme_c: c_float,
        vner_c: c_float,
        hmf2_c: c_float,
        hmf1_c: c_float,
        hme_c: c_float,
        hv1r_c: c_float,
        hv2r_c: c_float,
        hhalf_c: c_float,
        hxl_c: *mut c_float,
        scl_c: *mut c_float,
        amp_c: *mut c_float,
        iqu_c: *mut c_int,
    );

    pub fn xen_c(
        h_c: c_float,
        hmf2_c: c_float,
        nmf2_c: c_float,
        hme_c: c_float,
        n_c: c_int,
        hxl_c: *const c_float,
        scl_c: *const c_float,
        amp_c: *const c_float,
    ) -> c_float;

    pub fn teba_c(
        diplat_c: c_float,
        hour_c: c_float,
        season_c: c_int,
        tea_c: *mut c_float,
    );

    pub fn tede_c(
        h_c: c_float,
        xn_c: c_float,
        covsat_c: c_float,
    ) -> c_float;

    pub fn calion_c(
        invdip_c: c_float,
        mlt_c: c_float,
        alt_c: c_float,
        ddd_c: c_float,
        pf107_c: c_float,
        no_c: *mut c_float,
        nh_c: *mut c_float,
        nhe_c: *mut c_float,
        nn_c: *mut c_float,
    );

    pub fn chemion_c(
        jprint_c: c_int,
        alt_c: c_float,
        f107d_c: c_float,
        f10781_c: c_float,
        te_c: c_float,
        ti_c: c_float,
        tn_c: c_float,
        oxn_c: c_float,
        o2n_c: c_float,
        n2n_c: c_float,
        hen_c: c_float,
        hn_c: c_float,
        user_no_c: c_float,
        n4s_c: c_float,
        edens_c: c_float,
        user_oplus_c: c_float,
        szad_c: c_float,
        ro_c: *mut c_float,
        ro2_c: *mut c_float,
        rno_c: *mut c_float,
        rn2_c: *mut c_float,
        rn_c: *mut c_float,
        den_no_c: *mut c_float,
        den_n2d_c: *mut c_float,
        inewt_c: *mut c_int,
    );

    pub fn spreadf_brazil_c(
        daynr_c: c_int,
        idayy_c: c_int,
        f107d_c: c_float,
        lati_c: c_float,
        osfbr_c: *mut c_float,
    );

    pub fn clcmlt_c(
        iyyyy_c: c_int,
        ddd_c: c_int,
        uthr_c: c_float,
        glat_c: c_float,
        glon_c: c_float,
        mlt_c: *mut c_float,
    );

    pub fn b0_98_c(
        hour_c: c_float,
        sax_c: c_float,
        sux_c: c_float,
        nseasn_c: c_int,
        r_c: c_float,
        zlo_c: c_float,
        zmodip_c: c_float,
    ) -> c_float;
}

