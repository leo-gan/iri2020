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

end module iri_c_bindings
