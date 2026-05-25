module iri_c_bindings
    use iso_c_binding, only: c_int, c_float, c_bool, c_double

    implicit none

contains

    subroutine c_read_ig_rz() bind(C, name="read_ig_rz_c")
        call read_ig_rz()
    end subroutine c_read_ig_rz

    subroutine c_readapf107() bind(C, name="readapf107_c")
        call readapf107()
    end subroutine c_readapf107

    ! Signature for IRI_SUB
    ! SUBROUTINE IRI_SUB(JF,JMAG,ALATI,ALONG,IYYYY,MMDD,DHOUR,
    !    &    HEIBEG,HEIEND,HEISTP,OUTF,OARR)
    subroutine c_iri_sub(jf_c, jmag_c, alati_c, along_c, iyyyy_c, mmdd_c, dhour_c, &
                         heibeg_c, heiend_c, heistp_c, outf_c, oarr_c) bind(C, name="iri_sub_c")
        logical(c_bool), intent(in) :: jf_c(50)
        integer(c_int), intent(in), value :: jmag_c
        real(c_float), intent(in), value :: alati_c
        real(c_float), intent(in), value :: along_c
        integer(c_int), intent(in), value :: iyyyy_c
        integer(c_int), intent(in), value :: mmdd_c
        real(c_float), intent(in), value :: dhour_c
        real(c_float), intent(in), value :: heibeg_c
        real(c_float), intent(in), value :: heiend_c
        real(c_float), intent(in), value :: heistp_c
        real(c_float), intent(out) :: outf_c(20,1000)
        real(c_float), intent(out) :: oarr_c(100)

        logical :: jf(50)
        integer :: i

        do i = 1, 50
            jf(i) = jf_c(i)
        end do

        call IRI_SUB(jf, jmag_c, alati_c, along_c, iyyyy_c, mmdd_c, dhour_c, &
                     heibeg_c, heiend_c, heistp_c, outf_c, oarr_c)

    end subroutine c_iri_sub

    ! Signature for IRITEC
    ! SUBROUTINE IRITEC(ALATI,ALONG,JMAG,JF,IYYYY,MMDD,HOUR,
    !   &   HBGN,HEND,HSTEP,OARR,TEC,TECT)
    subroutine c_iritec(alati_c, along_c, jmag_c, jf_c, iyyyy_c, mmdd_c, hour_c, &
                        hbgn_c, hend_c, hstep_c, oarr_c, tec_c, tect_c) bind(C, name="iritec_c")
        real(c_float), intent(in), value :: alati_c
        real(c_float), intent(in), value :: along_c
        integer(c_int), intent(in), value :: jmag_c
        logical(c_bool), intent(in) :: jf_c(50)
        integer(c_int), intent(in), value :: iyyyy_c
        integer(c_int), intent(in), value :: mmdd_c
        real(c_float), intent(in), value :: hour_c
        real(c_float), intent(in), value :: hbgn_c
        real(c_float), intent(in), value :: hend_c
        real(c_float), intent(in), value :: hstep_c
        real(c_float), intent(inout) :: oarr_c(100)
        real(c_float), intent(out) :: tec_c
        real(c_float), intent(out) :: tect_c

        logical :: jf(50)
        integer :: i

        do i = 1, 50
            jf(i) = jf_c(i)
        end do

        call IRITEC(alati_c, along_c, jmag_c, jf, iyyyy_c, mmdd_c, hour_c, &
                    hbgn_c, hend_c, hstep_c, oarr_c, tec_c, tect_c)

    end subroutine c_iritec

    subroutine c_feldcof(year_c) bind(C, name="feldcof_c")
        real(c_float), intent(in), value :: year_c
        call FELDCOF(year_c)
    end subroutine c_feldcof

    subroutine c_feldg(glat_c, glon_c, alt_c, bnorth_c, beast_c, bdown_c, babs_c) bind(C, name="feldg_c")
        real(c_float), intent(in), value :: glat_c
        real(c_float), intent(in), value :: glon_c
        real(c_float), intent(in), value :: alt_c
        real(c_float), intent(out) :: bnorth_c
        real(c_float), intent(out) :: beast_c
        real(c_float), intent(out) :: bdown_c
        real(c_float), intent(out) :: babs_c
        call FELDG(glat_c, glon_c, alt_c, bnorth_c, beast_c, bdown_c, babs_c)
    end subroutine c_feldg

    subroutine c_igrf(iy_c, nm_c, r_c, t_c, f_c, br_c, bt_c, bf_c) bind(C, name="igrf_c")
        integer(c_int), intent(in), value :: iy_c
        integer(c_int), intent(in), value :: nm_c
        real(c_float), intent(in), value :: r_c
        real(c_float), intent(in), value :: t_c
        real(c_float), intent(in), value :: f_c
        real(c_float), intent(out) :: br_c
        real(c_float), intent(out) :: bt_c
        real(c_float), intent(out) :: bf_c
        call IGRF(iy_c, nm_c, r_c, t_c, f_c, br_c, bt_c, bf_c)
    end subroutine c_igrf

    subroutine c_igrf_dip(xlat_c, xlong_c, year_c, height_c, dec_c, dip_c, dipl_c, ymodip_c) bind(C, name="igrf_dip_c")
        real(c_float), intent(in), value :: xlat_c
        real(c_float), intent(in), value :: xlong_c
        real(c_float), intent(in), value :: year_c
        real(c_float), intent(in), value :: height_c
        real(c_float), intent(out) :: dec_c
        real(c_float), intent(out) :: dip_c
        real(c_float), intent(out) :: dipl_c
        real(c_float), intent(out) :: ymodip_c
        call igrf_dip(xlat_c, xlong_c, year_c, height_c, dec_c, dip_c, dipl_c, ymodip_c)
    end subroutine c_igrf_dip

    subroutine c_init_igrf() bind(C, name="init_igrf_c")
        real(c_float) :: era, aquad, bquad, dimo
        real(c_float) :: umr, pi
        real(c_float) :: argmax
        common /igrf1/ era, aquad, bquad, dimo
        common /const/ umr, pi
        common /argexp/ argmax
        
        era = 6371.2
        aquad = 6378.16 * 6378.16
        bquad = 6356.775 * 6356.775
        dimo = 0.311653
        pi = 3.141592653589793
        umr = pi / 180.0
        argmax = 87.3
    end subroutine c_init_igrf

    subroutine c_gtd7(iyd_c, sec_c, alt_c, glat_c, glong_c, stl_c, f107a_c, f107_c, ap_c, mass_c, d_c, t_c) bind(C, name="gtd7_c")
        integer(c_int), intent(in), value :: iyd_c
        real(c_float), intent(in), value :: sec_c
        real(c_float), intent(in), value :: alt_c
        real(c_float), intent(in), value :: glat_c
        real(c_float), intent(in), value :: glong_c
        real(c_float), intent(in), value :: stl_c
        real(c_float), intent(in), value :: f107a_c
        real(c_float), intent(in), value :: f107_c
        real(c_float), intent(in) :: ap_c(7)
        integer(c_int), intent(in), value :: mass_c
        real(c_float), intent(out) :: d_c(9)
        real(c_float), intent(out) :: t_c(2)

        real :: tlb, s, db04, db16, db28, db32, db40, db48, db01, za, t0, z0, g0, rl, dd, db14, tr12
        common/gts3c/tlb,s,db04,db16,db28,db32,db40,db48,db01,za,t0,z0,g0,rl,dd,db14,tr12
        real :: tn1(5), tn2(4), tn3(5), tgn1(2), tgn2(2), tgn3(2)
        common/meso7/tn1,tn2,tn3,tgn1,tgn2,tgn3

        call GTD7(iyd_c, sec_c, alt_c, glat_c, glong_c, stl_c, f107a_c, f107_c, ap_c, mass_c, d_c, t_c)

    end subroutine c_gtd7

    subroutine c_gtd7d(iyd_c, sec_c, alt_c, glat_c, glong_c, stl_c, f107a_c, f107_c, ap_c, mass_c, d_c, t_c) bind(C, name="gtd7d_c")
        integer(c_int), intent(in), value :: iyd_c
        real(c_float), intent(in), value :: sec_c
        real(c_float), intent(in), value :: alt_c
        real(c_float), intent(in), value :: glat_c
        real(c_float), intent(in), value :: glong_c
        real(c_float), intent(in), value :: stl_c
        real(c_float), intent(in), value :: f107a_c
        real(c_float), intent(in), value :: f107_c
        real(c_float), intent(in) :: ap_c(7)
        integer(c_int), intent(in), value :: mass_c
        real(c_float), intent(out) :: d_c(9)
        real(c_float), intent(out) :: t_c(2)
        call GTD7D(iyd_c, sec_c, alt_c, glat_c, glong_c, stl_c, f107a_c, f107_c, ap_c, mass_c, d_c, t_c)
    end subroutine c_gtd7d

    subroutine c_tselec(sv_c) bind(C, name="tselec_c")
        real(c_float), intent(in) :: sv_c(25)
        call TSELEC(sv_c)
    end subroutine c_tselec

    subroutine c_meters(meter_c) bind(C, name="meters_c")
        logical(c_bool), intent(in), value :: meter_c
        logical :: meter_val
        meter_val = meter_c
        call METERS(meter_val)
    end subroutine c_meters

    subroutine c_vfjmodelrocstart(vzm_c) bind(C, name="vfjmodelrocstart_c")
        real(c_float), intent(out) :: vzm_c(59, 25, 4, 11)
        call vfjmodelrocstart(vzm_c)
    end subroutine c_vfjmodelrocstart

    subroutine c_vfjmodelrocinit(f107_c, idoy_c, jseas_c, jsfl_c) bind(C, name="vfjmodelrocinit_c")
        real(c_float), intent(in), value :: f107_c
        integer(c_int), intent(in), value :: idoy_c
        integer(c_int), intent(out) :: jseas_c
        integer(c_int), intent(out) :: jsfl_c
        call vfjmodelrocinit(f107_c, idoy_c, jseas_c, jsfl_c)
    end subroutine c_vfjmodelrocinit

    subroutine c_vfjmodelroc(fjm_c, ttl_c, gglon_c, jseas_c, jsfl_c, viv_c) bind(C, name="vfjmodelroc_c")
        real(c_float), intent(in) :: fjm_c(59, 25, 4, 11)
        real(c_float), intent(in), value :: ttl_c
        real(c_float), intent(in), value :: gglon_c
        integer(c_int), intent(in), value :: jseas_c
        integer(c_int), intent(in), value :: jsfl_c
        real(c_float), intent(out) :: viv_c
        call vfjmodelroc(fjm_c, ttl_c, gglon_c, jseas_c, jsfl_c, viv_c)
    end subroutine c_vfjmodelroc
    subroutine c_get_igrz(aig_c, arz_c, iymst_c, iymend_c) bind(C, name="get_igrz_c")
        real(c_float), intent(out) :: aig_c(806)
        real(c_float), intent(out) :: arz_c(806)
        integer(c_int), intent(out) :: iymst_c
        integer(c_int), intent(out) :: iymend_c
        
        real :: aig(806), arz(806)
        integer :: iymst, iymend
        common /igrz/aig,arz,iymst,iymend
        
        aig_c = aig
        arz_c = arz
        iymst_c = iymst
        iymend_c = iymend
    end subroutine c_get_igrz

    subroutine c_get_apfa(aap_c, af107_c, n_c) bind(C, name="get_apfa_c")
        integer(c_int), intent(out) :: aap_c(27000, 9)
        real(c_float), intent(out) :: af107_c(27000, 3)
        integer(c_int), intent(out) :: n_c
        
        integer :: aap(27000,9)
        real :: af107(27000,3)
        integer :: n
        common /apfa/aap,af107,n
        
        aap_c = aap
        af107_c = af107
        n_c = n
    end subroutine c_get_apfa

    subroutine c_read_data_sd(month_c, coeff_month_c) bind(C, name="read_data_sd_c")
        integer(c_int), intent(in), value :: month_c
        real(c_double), intent(out) :: coeff_month_c(149, 48)
        
        double precision :: coeff_month(0:148, 0:47)
        
        call read_data_SD(month_c, coeff_month)
        
        coeff_month_c = coeff_month
    end subroutine c_read_data_sd

    subroutine c_read_coeff(month_c, is_ccir_c, f2_c, fm3_c) bind(C, name="read_coeff_c")
        integer(c_int), intent(in), value :: month_c
        logical(c_bool), intent(in), value :: is_ccir_c
        real(c_float), intent(out) :: f2_c(13, 76, 2)
        real(c_float), intent(out) :: fm3_c(9, 49, 2)
        
        character(256) filnam, prefix
        character(2) month_str
        integer :: iuccir
        
        iuccir = 10
        call get_data_prefix(prefix)
        write(month_str, '(I2)') month_c+10
        
        if (is_ccir_c) then
            filnam = trim(prefix) // 'ccir' // month_str // '.asc'
            open(iuccir, file=trim(filnam), status='old')
            read(iuccir, 4689) f2_c, fm3_c
            close(iuccir)
        else
            filnam = trim(prefix) // 'ursi' // month_str // '.asc'
            open(iuccir, file=trim(filnam), status='old')
            read(iuccir, 4689) f2_c
            close(iuccir)
        end if
        
4689    format(1X, 4E15.8)
    end subroutine c_read_coeff

    subroutine c_shellg(glat_c, glon_c, alt_c, fl_c, icode_c, b0_c) bind(C, name="shellg_c")
        real(c_float), intent(in), value :: glat_c
        real(c_float), intent(in), value :: glon_c
        real(c_float), intent(in), value :: alt_c
        real(c_float), intent(out) :: fl_c
        integer(c_int), intent(out) :: icode_c
        real(c_float), intent(out) :: b0_c
        call SHELLG(glat_c, glon_c, alt_c, fl_c, icode_c, b0_c)
    end subroutine c_shellg

    subroutine c_iondani(id_c, ismo_c, hx_c, zd_c, fd_c, fs_c, dion_c) bind(C, name="iondani_c")
        integer(c_int), intent(in), value :: id_c
        integer(c_int), intent(in), value :: ismo_c
        real(c_float), intent(in), value :: hx_c
        real(c_float), intent(in), value :: zd_c
        real(c_float), intent(in), value :: fd_c
        real(c_float), intent(in), value :: fs_c
        real(c_float), intent(out) :: dion_c(7)
        call iondani(id_c, ismo_c, hx_c, zd_c, fd_c, fs_c, dion_c)
    end subroutine c_iondani

    subroutine c_f1_prob(sza_c, glat_c, rz12_c, f1prob_c, f1probl_c) bind(C, name="f1_prob_c")
        real(c_float), intent(in), value :: sza_c
        real(c_float), intent(in), value :: glat_c
        real(c_float), intent(in), value :: rz12_c
        real(c_float), intent(out) :: f1prob_c
        real(c_float), intent(out) :: f1probl_c
        call f1_prob(sza_c, glat_c, rz12_c, f1prob_c, f1probl_c)
    end subroutine c_f1_prob

    subroutine c_set_xe_blocks(hmf2_c, xnmf2_c, hmf1_c, f1reg_c, &
                               b0_c, b1_c, c1_c, &
                               hz_c, t_c, hst_c, &
                               hme_c, xnme_c, hef_c, &
                               night_c, e_c, &
                               hmd_c, xnmd_c, hdx_c, &
                               d1_c, xkk_c, fp30_c, fp3u_c, fp1_c, fp2_c, &
                               beta_c, eta_c, delta_c, zeta_c, &
                               b2top_c, itopn_c, tcor1_c, tcor2_c) bind(C, name="set_xe_blocks_c")
        real(c_float), intent(in), value :: hmf2_c, xnmf2_c, hmf1_c
        logical(c_bool), intent(in), value :: f1reg_c
        real(c_float), intent(in), value :: b0_c, b1_c, c1_c
        real(c_float), intent(in), value :: hz_c, t_c, hst_c
        real(c_float), intent(in), value :: hme_c, xnme_c, hef_c
        logical(c_bool), intent(in), value :: night_c
        real(c_float), intent(in) :: e_c(4)
        real(c_float), intent(in), value :: hmd_c, xnmd_c, hdx_c
        real(c_float), intent(in), value :: d1_c, xkk_c, fp30_c, fp3u_c, fp1_c, fp2_c
        real(c_float), intent(in), value :: beta_c, eta_c, delta_c, zeta_c
        real(c_float), intent(in), value :: b2top_c
        integer(c_int), intent(in), value :: itopn_c
        real(c_float), intent(in), value :: tcor1_c, tcor2_c

        real :: HMF2, XNMF2, HMF1
        logical :: F1REG
        common /BLOCK1/ HMF2, XNMF2, HMF1, F1REG

        real :: B0, B1, C1
        common /BLOCK2/ B0, B1, C1

        real :: HZ, T_val, HST
        common /BLOCK3/ HZ, T_val, HST

        real :: HME, XNME, HEF
        common /BLOCK4/ HME, XNME, HEF

        logical :: NIGHT
        real :: E(4)
        common /BLOCK5/ NIGHT, E

        real :: HMD, XNMD, HDX
        common /BLOCK6/ HMD, XNMD, HDX

        real :: D1, XKK, FP30, FP3U, FP1, FP2
        common /BLOCK7/ D1, XKK, FP30, FP3U, FP1, FP2

        real :: BETA, ETA, DELTA, ZETA
        common /BLO10/ BETA, ETA, DELTA, ZETA

        real :: B2TOP, tcor1, tcor2
        integer :: itopn
        common /BLO11/ B2TOP, itopn, tcor1, tcor2

        HMF2 = hmf2_c
        XNMF2 = xnmf2_c
        HMF1 = hmf1_c
        F1REG = f1reg_c
        B0 = b0_c
        B1 = b1_c
        C1 = c1_c
        HZ = hz_c
        T_val = t_c
        HST = hst_c
        HME = hme_c
        XNME = xnme_c
        HEF = hef_c
        NIGHT = night_c
        E = e_c
        HMD = hmd_c
        XNMD = xnmd_c
        HDX = hdx_c
        D1 = d1_c
        XKK = xkk_c
        FP30 = fp30_c
        FP3U = fp3u_c
        FP1 = fp1_c
        FP2 = fp2_c
        BETA = beta_c
        ETA = eta_c
        DELTA = delta_c
        ZETA = zeta_c
        B2TOP = b2top_c
        itopn = itopn_c
        tcor1 = tcor1_c
        tcor2 = tcor2_c
    end subroutine c_set_xe_blocks

    subroutine c_get_xe_blocks(hmf2_c, xnmf2_c, hmf1_c, f1reg_c, &
                               b0_c, b1_c, c1_c, &
                               hz_c, t_c, hst_c, &
                               hme_c, xnme_c, hef_c, &
                               night_c, e_c, &
                               hmd_c, xnmd_c, hdx_c, &
                               d1_c, xkk_c, fp30_c, fp3u_c, fp1_c, fp2_c, &
                               beta_c, eta_c, delta_c, zeta_c, &
                               b2top_c, itopn_c, tcor1_c, tcor2_c) bind(C, name="get_xe_blocks_c")
        real(c_float), intent(out) :: hmf2_c, xnmf2_c, hmf1_c
        logical(c_bool), intent(out) :: f1reg_c
        real(c_float), intent(out) :: b0_c, b1_c, c1_c
        real(c_float), intent(out) :: hz_c, t_c, hst_c
        real(c_float), intent(out) :: hme_c, xnme_c, hef_c
        logical(c_bool), intent(out) :: night_c
        real(c_float), intent(out) :: e_c(4)
        real(c_float), intent(out) :: hmd_c, xnmd_c, hdx_c
        real(c_float), intent(out) :: d1_c, xkk_c, fp30_c, fp3u_c, fp1_c, fp2_c
        real(c_float), intent(out) :: beta_c, eta_c, delta_c, zeta_c
        real(c_float), intent(out) :: b2top_c
        integer(c_int), intent(out) :: itopn_c
        real(c_float), intent(out) :: tcor1_c, tcor2_c

        real :: HMF2, XNMF2, HMF1
        logical :: F1REG
        common /BLOCK1/ HMF2, XNMF2, HMF1, F1REG

        real :: B0, B1, C1
        common /BLOCK2/ B0, B1, C1

        real :: HZ, T_val, HST
        common /BLOCK3/ HZ, T_val, HST

        real :: HME, XNME, HEF
        common /BLOCK4/ HME, XNME, HEF

        logical :: NIGHT
        real :: E(4)
        common /BLOCK5/ NIGHT, E

        real :: HMD, XNMD, HDX
        common /BLOCK6/ HMD, XNMD, HDX

        real :: D1, XKK, FP30, FP3U, FP1, FP2
        common /BLOCK7/ D1, XKK, FP30, FP3U, FP1, FP2

        real :: BETA, ETA, DELTA, ZETA
        common /BLO10/ BETA, ETA, DELTA, ZETA

        real :: B2TOP, tcor1, tcor2
        integer :: itopn
        common /BLO11/ B2TOP, itopn, tcor1, tcor2

        hmf2_c = HMF2
        xnmf2_c = XNMF2
        hmf1_c = HMF1
        f1reg_c = F1REG
        b0_c = B0
        b1_c = B1
        c1_c = C1
        hz_c = HZ
        t_c = T_val
        hst_c = HST
        hme_c = HME
        xnme_c = XNME
        hef_c = HEF
        night_c = NIGHT
        e_c = E
        hmd_c = HMD
        xnmd_c = XNMD
        hdx_c = HDX
        d1_c = D1
        xkk_c = XKK
        fp30_c = FP30
        fp3u_c = FP3U
        fp1_c = FP1
        fp2_c = FP2
        beta_c = BETA
        eta_c = ETA
        delta_c = DELTA
        zeta_c = ZETA
        b2top_c = B2TOP
        itopn_c = itopn
        tcor1_c = tcor1
        tcor2_c = tcor2
    end subroutine c_get_xe_blocks

    real(c_float) function xe_1_c(h_c) bind(C, name="xe_1_c")
        real(c_float), intent(in), value :: h_c
        real :: XE_1
        xe_1_c = XE_1(h_c)
    end function xe_1_c

    real(c_float) function xe1_c(h_c) bind(C, name="xe1_c")
        real(c_float), intent(in), value :: h_c
        real :: XE1
        xe1_c = XE1(h_c)
    end function xe1_c

    real(c_float) function xe2_c(h_c) bind(C, name="xe2_c")
        real(c_float), intent(in), value :: h_c
        real :: XE2
        xe2_c = XE2(h_c)
    end function xe2_c

    real(c_float) function xe3_1_c(h_c) bind(C, name="xe3_1_c")
        real(c_float), intent(in), value :: h_c
        real :: XE3_1
        xe3_1_c = XE3_1(h_c)
    end function xe3_1_c

    real(c_float) function xe4_1_c(h_c) bind(C, name="xe4_1_c")
        real(c_float), intent(in), value :: h_c
        real :: XE4_1
        xe4_1_c = XE4_1(h_c)
    end function xe4_1_c

    real(c_float) function xe5_c(h_c) bind(C, name="xe5_c")
        real(c_float), intent(in), value :: h_c
        real :: XE5
        xe5_c = XE5(h_c)
    end function xe5_c

    real(c_float) function xe6_c(h_c) bind(C, name="xe6_c")
        real(c_float), intent(in), value :: h_c
        real :: XE6
        xe6_c = XE6(h_c)
    end function xe6_c

    real(c_float) function dxe1n_c(h_c) bind(C, name="dxe1n_c")
        real(c_float), intent(in), value :: h_c
        real :: DXE1N
        dxe1n_c = DXE1N(h_c)
    end function dxe1n_c

    real(c_float) function topq_c(h_c, no_c, hmax_c, ho_c) bind(C, name="topq_c")
        real(c_float), intent(in), value :: h_c
        real(c_float), intent(in), value :: no_c
        real(c_float), intent(in), value :: hmax_c
        real(c_float), intent(in), value :: ho_c
        real :: TOPQ
        topq_c = TOPQ(h_c, no_c, hmax_c, ho_c)
    end function topq_c

    real(c_float) function zero_c(delta_c) bind(C, name="zero_c")
        real(c_float), intent(in), value :: delta_c
        real :: ZERO
        zero_c = ZERO(delta_c)
    end function zero_c

    real(c_float) function foeedi_c(cov_c, xhi_c, xhim_c, xlati_c) bind(C, name="foeedi_c")
        real(c_float), intent(in), value :: cov_c
        real(c_float), intent(in), value :: xhi_c
        real(c_float), intent(in), value :: xhim_c
        real(c_float), intent(in), value :: xlati_c
        real :: FOEEDI
        foeedi_c = FOEEDI(cov_c, xhi_c, xhim_c, xlati_c)
    end function foeedi_c

    real(c_float) function xmded_c(xhi_c, r_c, yw_c) bind(C, name="xmded_c")
        real(c_float), intent(in), value :: xhi_c
        real(c_float), intent(in), value :: r_c
        real(c_float), intent(in), value :: yw_c
        real :: XMDED
        xmded_c = XMDED(xhi_c, r_c, yw_c)
    end function xmded_c

    subroutine valgul_c(xhi_c, hvb_c, vwu_c, vwa_c, vdp_c) bind(C, name="valgul_c")
        real(c_float), intent(in), value :: xhi_c
        real(c_float), intent(out) :: hvb_c, vwu_c, vwa_c, vdp_c
        call VALGUL(xhi_c, hvb_c, vwu_c, vwa_c, vdp_c)
    end subroutine valgul_c

    subroutine spharm_c(c_c, l_c, m_c, colat_c, az_c) bind(C, name="spharm_c")
        real(c_float), intent(inout) :: c_c(82)
        integer(c_int), intent(in), value :: l_c, m_c
        real(c_float), intent(in), value :: colat_c, az_c
        call SPHARM(c_c, l_c, m_c, colat_c, az_c)
    end subroutine spharm_c

    subroutine spharm_ik_c(c_c, l_c, m_c, colat_c, az_c) bind(C, name="spharm_ik_c")
        real(c_float), intent(inout) :: c_c(82)
        integer(c_int), intent(in), value :: l_c, m_c
        real(c_float), intent(in), value :: colat_c, az_c
        call SPHARM_IK(c_c, l_c, m_c, colat_c, az_c)
    end subroutine spharm_ik_c

    subroutine dregion_c(z_c, it_c, f_c, vkp_c, f5sw_c, f6wa_c, elg_c) bind(C, name="dregion_c")
        real(c_float), intent(in), value :: z_c
        integer(c_int), intent(in), value :: it_c
        real(c_float), intent(in), value :: f_c, vkp_c, f5sw_c, f6wa_c
        real(c_float), intent(out) :: elg_c(7)
        call DRegion(z_c, it_c, f_c, vkp_c, f5sw_c, f6wa_c, elg_c)
    end subroutine dregion_c

    subroutine f00_c(hgt_c, glat1_c, iday_c, zang_c, f107t_c, edens_c, ierror_c) bind(C, name="f00_c")
        real(c_float), intent(in), value :: hgt_c, glat1_c
        integer(c_int), intent(in), value :: iday_c
        real(c_float), intent(in), value :: zang_c, f107t_c
        real(c_float), intent(out) :: edens_c
        integer(c_int), intent(out) :: ierror_c
        call F00(hgt_c, glat1_c, iday_c, zang_c, f107t_c, edens_c, ierror_c)
    end subroutine f00_c

    real(c_float) function shamdb0d_c(rlat_c, flon_c, t_c, rz_c) bind(C, name="shamdb0d_c")
        real(c_float), intent(in), value :: rlat_c, flon_c, t_c, rz_c
        real :: b
        call SHAMDB0D(rlat_c, flon_c, t_c, rz_c, b)
        shamdb0d_c = b
    end function shamdb0d_c

    real(c_float) function shab1d_c(flat_c, flon_c, t_c, rz_c) bind(C, name="shab1d_c")
        real(c_float), intent(in), value :: flat_c, flon_c, t_c, rz_c
        real :: b
        call SHAB1D(flat_c, flon_c, t_c, rz_c, b)
        shab1d_c = b
    end function shab1d_c

end module iri_c_bindings


