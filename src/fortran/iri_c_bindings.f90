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
        common /igrf1/ era, aquad, bquad, dimo
        common /const/ umr, pi
        
        era = 6371.2
        aquad = 6378.16 * 6378.16
        bquad = 6356.775 * 6356.775
        dimo = 0.311653
        pi = 3.141592653589793
        umr = pi / 180.0
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

end module iri_c_bindings
