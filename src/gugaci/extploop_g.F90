!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!***********************************************************************
! Module_calcu  Completed LOOP multiply MOs and partly complete DM1 and

! SUBROUTINE lp10_arbrbr_ext_calcuvalue_G
! SUBROUTINE lp11_arblbr_ext_calcuvalue_G
! SUBROUTINE lp12_arblbl_ext_calcuvalue_G
! SUBROUTINE lp9_drlbl_ext_calcuvalue_wyb_G
! SUBROUTINE lp_ar_coe_sd_calcuvalue_wyb_G
! SUBROUTINE lp_ar_coe_td_calcuvalue_wyb_G
! SUBROUTINE lp_ar_coe_dv_calcuvalue_wyb_G
! SUBROUTINE lp_arbr_ext_svtv_calcuvalue_wyb_G
! SUBROUTINE lp_drr_ext_svtv_calcuvalue_wyb_G
! SUBROUTINE lp_arbl_ext_ss_calcuvalue_G
! SUBROUTINE lp_arbl_ext_st_calcuvalue_G
! SUBROUTINE lp_arbl_ext_ts_calcuvalue_G
! SUBROUTINE lp_arbl_ext_tt_calcuvalue_G
! SUBROUTINE lp_arbl_ext_dd_calcuvalue_G
! SUBROUTINE lp_drl_ext_SS_calcuvalue_G
! SUBROUTINE lp_drl_ext_ST_calcuvalue_G
! SUBROUTINE lp_drl_ext_TS_calcuvalue_G
! SUBROUTINE lp_drl_ext_TT_calcuvalue_G
! SUBROUTINE lp_drl_ext_dd_calcuvalue_wyb_G
! SUBROUTINE lp678_ext_wyb_calcuvalue_G_1
! SUBROUTINE lp678_ext_wyb_calcuvalue_G_2
! SUBROUTINE lp_drl_sum_SS_calcuvalue_G
! SUBROUTINE lp_drl_sum_TT_calcuvalue_G
! SUBROUTINE lp8_drlbr_sum_calcuvalue_wyb_G
! SUBROUTINE lp9_drlbl_sum_calcuvalue_wyb_G
! SUBROUTINE gsd_ext_sequence_G

subroutine lp_drl_ext_SS_calcuvalue_G(lri,nlp_value)

use gugaci_global, only: ibsm_ext, iesm_ext, index_lpext, index_lpext1, logic_g2g4a, ng_sm, norb_number, value_lpext, &
                         value_lpext1, w0_plp, w0g2a, w0g36a, w1_plp, w1g2a, w1g36a
use Constants, only: Two
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: lri
integer(kind=iwp), intent(out) :: nlp_value
integer(kind=iwp) :: ibend, ibsta, ira, irb, ivalue, lmb, lra, lrb, nxo
real(kind=wp) :: w0lp, w1lp

ivalue = 0

! G2G4b
if (logic_g2g4a) then
  w0lp = w0_plp*w0g2a
  w1lp = w1_plp*w1g2a

  do lmb=1,ng_sm
    ibsta = ibsm_ext(lmb)
    ibend = iesm_ext(lmb)
    do irb=ibsta,ibend
      lrb = norb_number(irb)
      do ira=ibsta,irb-1
        lra = norb_number(ira)

        ivalue = ivalue+2
        ! Drl -- B^lA^r =4_3
        call TRANS_IJKL_INTPOS(lra,lrb,lri,lri,NXO)
        index_lpext(ivalue) = NXO
        value_lpext(ivalue) = -Two*w0lp

        call TRANS_IJKL_INTPOS(lra,lri,lrb,lri,NXO)
        index_lpext1(ivalue) = NXO
        value_lpext1(ivalue) = (w0lp-w1lp)

      end do
    end do
  end do
end if

w0lp = w0_plp*w0g36a
w1lp = w1_plp*w1g36a

do lmb=1,ng_sm
  ibsta = ibsm_ext(lmb)
  ibend = iesm_ext(lmb)
  do irb=ibsta,ibend
    lrb = norb_number(irb)
    do ira=ibsta,irb-1
      lra = norb_number(ira)
      ivalue = ivalue+1

      call TRANS_IJKL_INTPOS(lra,lrb,lri,lri,NXO)
      index_lpext(ivalue) = NXO
      value_lpext(ivalue) = -Two*w0lp

      call TRANS_IJKL_INTPOS(lra,lri,lrb,lri,NXO)
      index_lpext1(ivalue) = NXO
      value_lpext1(ivalue) = (w0lp-w1lp)

    end do
  end do
end do
nlp_value = ivalue

end subroutine lp_drl_ext_SS_calcuvalue_G

subroutine lp_drl_ext_ST_calcuvalue_G(lri,nlp_value)

use gugaci_global, only: ibsm_ext, iesm_ext, index_lpext, index_lpext1, ism_g1415, logic_g1415, logic_g2g4b, ng_sm, norb_number, &
                         value_lpext, value_lpext1, w1_plp, w1g14a, w1g36a, w1g4b
use Symmetry_Info, only: Mul
use Constants, only: Two
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: lri
integer(kind=iwp), intent(out) :: nlp_value
integer(kind=iwp) :: ia, iaend, iasta, ib, ibend, ibsta, ira, irb, isma, ismb, ivalue, lmb, lra, lrb, nxo
real(kind=wp) :: w1lp

ivalue = 0
! G1415
if (logic_g1415) then

  w1lp = w1_plp*w1g14a
  !=========================lyb=========================================
  ! THE W1LP SHOULD BE MULTIPLED BY TWO,BECAUSE WE JUST USE HALF OF THEM
  ! IS H*C CALCULATIONS, ???

  w1lp = w1lp*Two

  do ismb=1,ng_sm
    isma = Mul(ismb,ism_g1415)
    if (isma > ismb) cycle
    ibsta = ibsm_ext(ismb)
    ibend = iesm_ext(ismb)
    iasta = ibsm_ext(isma)
    iaend = iesm_ext(isma)
    if (ismb == isma) ibsta = ibsta+1
    do ib=ibsta,ibend
      lrb = norb_number(ib)
      do ia=iasta,min(iaend,ib-1)
        lra = norb_number(ia)

        ivalue = ivalue+1
        call TRANS_IJKL_INTPOS(lra,LRI,lra,LRI,NXO)
        index_lpext(ivalue) = NXO
        value_lpext(ivalue) = w1lp

        call TRANS_IJKL_INTPOS(lrb,LRI,lrb,LRI,NXO)

        index_lpext1(ivalue) = NXO
        value_lpext1(ivalue) = -w1lp

      end do
    end do
  end do
end if

! G2G4b
if (logic_g2g4b) then
  w1lp = w1_plp*w1g4b
  do lmb=1,ng_sm
    ibsta = ibsm_ext(lmb)
    ibend = iesm_ext(lmb)
    do irb=ibsta,ibend
      lrb = norb_number(irb)
      do ira=ibsta,irb-1
        lra = norb_number(ira)

        ivalue = ivalue+1
        call TRANS_IJKL_INTPOS(lra,lri,lrb,lri,NXO)
        index_lpext(ivalue) = NXO
        value_lpext(ivalue) = -w1lp

        ivalue = ivalue+1
        index_lpext(ivalue) = NXO
        value_lpext(ivalue) = w1lp

      end do
    end do
  end do
end if

w1lp = w1_plp*w1g36a
do lmb=1,ng_sm
  ibsta = ibsm_ext(lmb)
  ibend = iesm_ext(lmb)
  do irb=ibsta,ibend
    lrb = norb_number(irb)
    do ira=ibsta,irb-1
      lra = norb_number(ira)
      ivalue = ivalue+1
      call TRANS_IJKL_INTPOS(lra,lri,lrb,lri,NXO)
      index_lpext(ivalue) = NXO
      value_lpext(ivalue) = -w1lp
    end do
  end do
end do
nlp_value = ivalue

end subroutine lp_drl_ext_ST_calcuvalue_G

subroutine lp_drl_ext_TT_calcuvalue_G(lri,n1415_value,nlp_value)

use gugaci_global, only: ibsm_ext, iesm_ext, index_lpext, index_lpext1, ism_g1415, logic_g1415, ng_sm, norb_number, value_lpext, &
                         value_lpext1, w0_plp, w0g14a, w0g15a, w0g36a, w1_plp, w1g14a, w1g15a, w1g36a
use Symmetry_Info, only: Mul
use Constants, only: Two
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: lri
integer(kind=iwp), intent(out) :: n1415_value, nlp_value
integer(kind=iwp) :: ia, iaend, iasta, ib, ibend, ibsta, ira, irb, isma, ismb, ivalue, lmb, lra, lrb, nxo
real(kind=wp) :: w014, w015, w0lp, w114, w115, w14lp, w15lp, w1lp

ivalue = 0
! G1415
if (logic_g1415) then

  w014 = w0_plp*w0g14a
  w114 = w1_plp*w1g14a
  w015 = w0_plp*w0g15a
  w115 = w1_plp*w1g15a

  w14lp = w014-w114
  w15lp = w015-w115
  !=========================lyb=========================================
  ! THE W1LP SHOULD BE MULTIPLED BY TWO,BECAUSE WE JUST USE HALF OF THEM
  ! IS H*C CALCULATIONS, ???

  w14lp = w14lp*Two
  w15lp = w15lp*Two

  do ismb=1,ng_sm
    isma = Mul(ismb,ism_g1415)
    if (isma > ismb) cycle
    ibsta = ibsm_ext(ismb)
    ibend = iesm_ext(ismb)
    iasta = ibsm_ext(isma)
    iaend = iesm_ext(isma)
    if (ismb == isma) ibsta = ibsta+1
    do ib=ibsta,ibend
      lrb = norb_number(ib)
      do ia=iasta,min(iaend,ib-1)
        lra = norb_number(ia)
        ivalue = ivalue+1
        call TRANS_IJKL_INTPOS(lra,LRI,lra,LRI,NXO)
        index_lpext(ivalue) = NXO
        value_lpext(ivalue) = w15lp
        call TRANS_IJKL_INTPOS(lrb,LRI,lrb,LRI,NXO)
        index_lpext1(ivalue) = NXO
        value_lpext1(ivalue) = w14lp
      end do
    end do
  end do
end if

n1415_value = ivalue

w0lp = w0_plp*w0g36a
w1lp = w1_plp*w1g36a

do lmb=1,ng_sm
  ibsta = ibsm_ext(lmb)
  ibend = iesm_ext(lmb)
  do irb=ibsta,ibend
    lrb = norb_number(irb)
    do ira=ibsta,irb-1
      lra = norb_number(ira)
      ivalue = ivalue+1

      call TRANS_IJKL_INTPOS(lra,lrb,lri,lri,NXO)
      index_lpext(ivalue) = NXO
      value_lpext(ivalue) = -Two*w0lp

      call TRANS_IJKL_INTPOS(lra,lri,lrb,lri,NXO)
      index_lpext1(ivalue) = NXO
      value_lpext1(ivalue) = (w0lp-w1lp)

    end do
  end do
end do
nlp_value = ivalue

end subroutine lp_drl_ext_TT_calcuvalue_G

subroutine lp_drl_SUM_TT_calcuvalue_G(lri,n1415_value,nlp_value)

use gugaci_global, only: ibsm_ext, iesm_ext, index_lpext, index_lpext1, ism_g1415, logic_g1415, ng_sm, norb_number, value_lpext, &
                         value_lpext1, w0_plp, w0g14a, w0g15a, w0g36a, w1_plp, w1g14a, w1g15a
use Symmetry_Info, only: Mul
use Constants, only: Two
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: lri
integer(kind=iwp), intent(out) :: n1415_value, nlp_value
integer(kind=iwp) :: ia, iaend, iasta, ib, ibend, ibsta, ira, irb, isma, ismb, ivalue, lmb, lra, lrb, nxo
real(kind=wp) :: w014, w015, w0lp, w114, w115, w14lp, w15lp

ivalue = 0
! G1415
if (logic_g1415) then

  w014 = w0_plp*w0g14a
  w114 = w1_plp*w1g14a
  w015 = w0_plp*w0g15a
  w115 = w1_plp*w1g15a

  w14lp = w014-w114
  w15lp = w015-w115
  !=========================lyb=========================================
  ! THE W1LP SHOULD BE MULTIPLED BY TWO,BECAUSE WE JUST USE HALF OF THEM
  ! IS H*C CALCULATIONS, ???

  w14lp = w14lp*Two
  w15lp = w15lp*Two
  do ismb=1,ng_sm
    isma = Mul(ismb,ism_g1415)
    if (isma > ismb) cycle
    ibsta = ibsm_ext(ismb)
    ibend = iesm_ext(ismb)
    iasta = ibsm_ext(isma)
    iaend = iesm_ext(isma)
    if (ismb == isma) ibsta = ibsta+1
    do ib=ibsta,ibend
      lrb = norb_number(ib)
      do ia=iasta,min(iaend,ib-1)
        lra = norb_number(ia)
        ivalue = ivalue+1
        call TRANS_IJKL_INTPOS(lra,LRI,lra,LRI,NXO)
        index_lpext(ivalue) = NXO
        value_lpext(ivalue) = w15lp
        call TRANS_IJKL_INTPOS(lrb,LRI,lrb,LRI,NXO)
        index_lpext1(ivalue) = NXO
        value_lpext1(ivalue) = w14lp
      end do
    end do
  end do
end if

n1415_value = ivalue

w0lp = w0_plp*w0g36a

do lmb=1,ng_sm
  ibsta = ibsm_ext(lmb)
  ibend = iesm_ext(lmb)
  do irb=ibsta,ibend
    lrb = norb_number(irb)
    do ira=ibsta,irb-1
      lra = norb_number(ira)
      ivalue = ivalue+1
      call TRANS_IJKL_INTPOS(lra,lrb,lri,lri,NXO)
      index_lpext(ivalue) = NXO
      value_lpext(ivalue) = -Two*w0lp
      call TRANS_IJKL_INTPOS(lra,lri,lrb,lri,NXO)
      index_lpext1(ivalue) = NXO
      value_lpext1(ivalue) = w0lp
    end do
  end do
end do
nlp_value = ivalue

end subroutine lp_drl_SUM_TT_calcuvalue_G

subroutine lp_drl_ext_TS_calcuvalue_G(lri,nlp_value)

use gugaci_global, only: ibsm_ext, iesm_ext, index_lpext, index_lpext1, ism_g1415, logic_g1415, logic_g2g4a, ng_sm, norb_number, &
                         value_lpext, value_lpext1, w1_plp, w1g14a, w1g2a, w1g36a
use Symmetry_Info, only: Mul
use Constants, only: Two
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: lri
integer(kind=iwp), intent(out) :: nlp_value
integer(kind=iwp) :: ia, iaend, iasta, ib, ibend, ibsta, ira, irb, isma, ismb, ivalue, lmb, lra, lrb, nxo
real(kind=wp) :: w1lp

ivalue = 0
! G1415
if (logic_g1415) then

  w1lp = w1_plp*w1g14a
  !=========================lyb=========================================
  ! THE W1LP SHOULD BE MULTIPLED BY TWO,BECAUSE WE JUST USE HALF OF THEM
  ! IS H*C CALCULATIONS, ???

  w1lp = w1lp*Two
  do ismb=1,ng_sm
    isma = Mul(ismb,ism_g1415)
    if (isma > ismb) cycle
    ibsta = ibsm_ext(ismb)
    ibend = iesm_ext(ismb)
    iasta = ibsm_ext(isma)
    iaend = iesm_ext(isma)
    if (ismb == isma) ibsta = ibsta+1
    do ib=ibsta,ibend
      lrb = norb_number(ib)
      do ia=iasta,min(iaend,ib-1)
        lra = norb_number(ia)
        ivalue = ivalue+1
        call TRANS_IJKL_INTPOS(lra,LRI,lra,LRI,NXO)
        index_lpext(ivalue) = NXO
        value_lpext(ivalue) = w1lp
        call TRANS_IJKL_INTPOS(lrb,LRI,lrb,LRI,NXO)
        index_lpext1(ivalue) = NXO
        value_lpext1(ivalue) = -w1lp
      end do
    end do
  end do
end if

! G2G4b
if (logic_g2g4a) then
  !w0lp = w0_plp*w0g2a
  w1lp = w1_plp*w1g2a

  !valuelptmp1 = w0lp
  !w0lp = w0lp-w1lp
  !w1lp = -valuelptmp1*Two
  !valuelptmp1 = ww0lp
  !ww0lp = ww0lp-ww1lp
  !ww1lp = -valuelptmp1*Two

  do lmb=1,ng_sm
    ibsta = ibsm_ext(lmb)
    ibend = iesm_ext(lmb)
    do irb=ibsta,ibend
      lrb = norb_number(irb)
      do ira=ibsta,irb-1
        lra = norb_number(ira)

        ivalue = ivalue+1
        call TRANS_IJKL_INTPOS(lra,lri,lrb,lri,NXO)
        index_lpext(ivalue) = NXO
        value_lpext(ivalue) = -w1lp

        ivalue = ivalue+1
        index_lpext(ivalue) = NXO
        value_lpext(ivalue) = w1lp
      end do
    end do
  end do
end if

w1lp = w1_plp*w1g36a

do lmb=1,ng_sm
  ibsta = ibsm_ext(lmb)
  ibend = iesm_ext(lmb)
  do irb=ibsta,ibend
    lrb = norb_number(irb)
    do ira=ibsta,irb-1
      lra = norb_number(ira)
      ivalue = ivalue+1
      call TRANS_IJKL_INTPOS(lra,lri,lrb,lri,NXO)
      index_lpext(ivalue) = NXO
      value_lpext(ivalue) = -w1lp
    end do
  end do
end do
nlp_value = ivalue

end subroutine lp_drl_ext_TS_calcuvalue_G

subroutine lp_arbl_ext_st_calcuvalue_G(lri,lrj,nlp_value)

use gugaci_global, only: ibsm_ext, iesm_ext, index_lpext, index_lpext1, ism_g1415, logic_g13, logic_g1415, logic_g2g4a, &
                         logic_g2g4b, lsm_inn, ng_sm, norb_ext, norb_number, value_lpext, value_lpext1, w0_plp, w0g13a, w0g14a, &
                         w0g15a, w0g2a, w0g2b, w0g36a, w0g36b, w0g4a, w0g4b, w1_plp, w1g14a, w1g15a, w1g2a, w1g2b, w1g36a, w1g36b, &
                         w1g4a, w1g4b
use Symmetry_Info, only: Mul
use Constants, only: Two
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: lri, lrj
integer(kind=iwp), intent(out) :: nlp_value
integer(kind=iwp) :: ia, iaend, iasta, ib, ibend, ibsta, isma, ismb, ivalue, lra, lrb, lsma, lsmb, lsmi, lsmij, lsmj, nxo
real(kind=wp) :: valuelptmp1, w0lp, w1lp, ww0lp, ww1lp

ivalue = 0
! G1415
if (logic_g1415) then

  w0lp = w0_plp*w0g14a
  w1lp = w1_plp*w1g14a
  ww0lp = w0_plp*w0g15a
  ww1lp = w1_plp*w1g15a

  valuelptmp1 = w0lp
  w0lp = w0lp-w1lp
  w1lp = -valuelptmp1*Two
  valuelptmp1 = ww0lp
  ww0lp = ww0lp-ww1lp
  ww1lp = -valuelptmp1*Two

  do ismb=1,ng_sm
    isma = Mul(ismb,ism_g1415)
    if (isma > ismb) cycle
    ibsta = ibsm_ext(ismb)
    ibend = iesm_ext(ismb)
    iasta = ibsm_ext(isma)
    iaend = iesm_ext(isma)
    if (ismb == isma) ibsta = ibsta+1
    do ib=ibsta,ibend
      lrb = norb_number(ib)
      do ia=iasta,min(iaend,ib-1)
        lra = norb_number(ia)

        ivalue = ivalue+1
        call TRANS_IJKL_INTPOS(LRJ,lrb,LRI,lrb,NXO)
        index_lpext(ivalue) = NXO
        value_lpext(ivalue) = w0lp
        call TRANS_IJKL_INTPOS(LRJ,LRI,lrb,lrb,NXO)
        index_lpext1(ivalue) = NXO
        value_lpext1(ivalue) = w1lp

        ivalue = ivalue+1
        call TRANS_IJKL_INTPOS(LRJ,lra,LRI,lra,NXO)
        index_lpext(ivalue) = NXO
        value_lpext(ivalue) = ww0lp
        call TRANS_IJKL_INTPOS(LRJ,LRI,lra,lra,NXO)
        index_lpext1(ivalue) = NXO
        value_lpext1(ivalue) = ww1lp
      end do
    end do
  end do
end if
! G13
if (logic_g13) then
  w0lp = w0g13a*w0_plp
  do ia=1,norb_ext
    lra = norb_number(ia)

    ivalue = ivalue+1
    call TRANS_IJKL_INTPOS(LRJ,lra,LRI,lra,NXO)
    index_lpext(ivalue) = NXO
    value_lpext(ivalue) = w0lp
    call TRANS_IJKL_INTPOS(LRJ,LRI,lra,lra,NXO)
    index_lpext1(ivalue) = NXO
    value_lpext1(ivalue) = -w0lp*Two

    ivalue = ivalue+1
    index_lpext(ivalue) = 0
    index_lpext1(ivalue) = 0

  end do
end if

lsmi = lsm_inn(lri)
lsmj = lsm_inn(lrj)
lsmij = Mul(lsmi,lsmj)

! G2G4a
if (logic_g2g4a) then

  w0lp = w0_plp*w0g2a
  w1lp = w1_plp*w1g2a
  ww0lp = w0_plp*w0g4a
  ww1lp = w1_plp*w1g4a

  valuelptmp1 = w0lp
  w0lp = w0lp-w1lp
  w1lp = -valuelptmp1*Two
  valuelptmp1 = ww0lp
  ww0lp = ww0lp-ww1lp
  ww1lp = -valuelptmp1*Two

  do lsmb=1,ng_sm
    lsma = Mul(lsmij,lsmb)
    if (lsma > lsmb) cycle
    ibsta = ibsm_ext(lsmb)
    ibend = iesm_ext(lsmb)
    iasta = ibsm_ext(lsma)
    iaend = iesm_ext(lsma)
    if (lsmb == lsma) ibsta = ibsta+1
    do ib=ibsta,ibend
      LRB = norb_number(ib)
      do ia=iasta,min(iaend,ib-1)
        LRA = norb_number(ia)
        ! ArBl -- B^lA^r =12
        ivalue = ivalue+1
        call TRANS_IJKL_INTPOS(lra,lri,lrj,lrb,NXO)
        index_lpext(ivalue) = NXO
        value_lpext(ivalue) = w0lp
        call TRANS_IJKL_INTPOS(lra,lrb,lrj,lri,NXO)
        index_lpext1(ivalue) = NXO
        value_lpext1(ivalue) = w1lp
        ! ArBl -- B^rA^l =11
        ivalue = ivalue+1
        call TRANS_IJKL_INTPOS(lra,lrj,lrb,lri,NXO)
        index_lpext(ivalue) = NXO
        value_lpext(ivalue) = ww0lp
        call TRANS_IJKL_INTPOS(lra,lrb,lrj,lri,NXO)
        index_lpext1(ivalue) = NXO
        value_lpext1(ivalue) = ww1lp
      end do
    end do
  end do
else

  ! G2G4b
  if (logic_g2g4b) then

    w0lp = w0_plp*w0g2b
    w1lp = w1_plp*w1g2b
    ww0lp = w0_plp*w0g4b
    ww1lp = w1_plp*w1g4b

    valuelptmp1 = w0lp
    w0lp = w0lp-w1lp
    w1lp = -valuelptmp1*Two
    valuelptmp1 = ww0lp
    ww0lp = ww0lp-ww1lp
    ww1lp = -valuelptmp1*Two

    do lsmb=1,ng_sm
      lsma = Mul(lsmij,lsmb)
      if (lsma > lsmb) cycle
      ibsta = ibsm_ext(lsmb)
      ibend = iesm_ext(lsmb)
      iasta = ibsm_ext(lsma)
      iaend = iesm_ext(lsma)
      if (lsmb == lsma) ibsta = ibsta+1
      do ib=ibsta,ibend
        LRB = norb_number(ib)
        do ia=iasta,min(iaend,ib-1)
          LRA = norb_number(ia)
          ! ArBl -- B^lA^r =12
          ivalue = ivalue+1
          call TRANS_IJKL_INTPOS(lra,lri,lrj,lrb,NXO)
          index_lpext(ivalue) = NXO
          value_lpext(ivalue) = ww0lp
          call TRANS_IJKL_INTPOS(lra,lrb,lrj,lri,NXO)
          index_lpext1(ivalue) = NXO
          value_lpext1(ivalue) = ww1lp
          ! ArBl -- B^rA^l =11
          ivalue = ivalue+1
          call TRANS_IJKL_INTPOS(lra,lrj,lrb,lri,NXO)
          index_lpext(ivalue) = NXO
          value_lpext(ivalue) = w0lp
          call TRANS_IJKL_INTPOS(lra,lrb,lrj,lri,NXO)
          index_lpext1(ivalue) = NXO
          value_lpext1(ivalue) = w1lp
        end do
      end do
    end do

  end if

end if
! G36a

w0lp = w0_plp*w0g36a
w1lp = w1_plp*w1g36a
valuelptmp1 = w0lp
w0lp = w0lp-w1lp
w1lp = -valuelptmp1*Two
do lsmb=1,ng_sm
  lsma = Mul(lsmij,lsmb)
  if (lsma > lsmb) cycle
  ibsta = ibsm_ext(lsmb)
  ibend = iesm_ext(lsmb)
  iasta = ibsm_ext(lsma)
  iaend = iesm_ext(lsma)
  if (lsmb == lsma) ibsta = ibsta+1
  do ib=ibsta,ibend
    LRB = norb_number(ib)
    do ia=iasta,min(iaend,ib-1)
      LRA = norb_number(ia)
      ! ArBl -- B^lA^r =12
      ivalue = ivalue+1
      call TRANS_IJKL_INTPOS(lra,lri,lrj,lrb,NXO)
      index_lpext(ivalue) = NXO
      value_lpext(ivalue) = w0lp
      call TRANS_IJKL_INTPOS(lra,lrb,lrj,lri,NXO)
      index_lpext1(ivalue) = NXO
      value_lpext1(ivalue) = w1lp
    end do
  end do
end do
! G36b

w0lp = w0_plp*w0g36b
w1lp = w1_plp*w1g36b
valuelptmp1 = w0lp
w0lp = w0lp-w1lp
w1lp = -valuelptmp1*Two
do lsmb=1,ng_sm
  lsma = Mul(lsmij,lsmb)
  if (lsma > lsmb) cycle
  ibsta = ibsm_ext(lsmb)
  ibend = iesm_ext(lsmb)
  iasta = ibsm_ext(lsma)
  iaend = iesm_ext(lsma)
  if (lsmb == lsma) ibsta = ibsta+1
  do ib=ibsta,ibend
    LRB = norb_number(ib)
    do ia=iasta,min(iaend,ib-1)
      LRA = norb_number(ia)
      ivalue = ivalue+1
      ! ArBl -- B^rA^l =11
      call TRANS_IJKL_INTPOS(lra,lrj,lrb,lri,NXO)
      index_lpext(ivalue) = NXO
      value_lpext(ivalue) = w0lp
      call TRANS_IJKL_INTPOS(lra,lrb,lrj,lri,NXO)
      index_lpext1(ivalue) = NXO
      value_lpext1(ivalue) = w1lp
    end do
  end do
end do
nlp_value = ivalue

end subroutine lp_arbl_ext_st_calcuvalue_G

subroutine lp10_arbrbr_ext_calcuvalue_G(intentry,isma,nlp_value)

use gugaci_global, only: ibsm_ext, index_lpext, index_lpext1, intind_ijka, lsm_inn, ngw2, ngw3, nlsm_ext, norb_frz, norb_inn, &
                         norb_number, value_lpext, value_lpext1, w0_sdplp, w0_sdplp25, w0g25, w1_sdplp, w1_sdplp25
use Symmetry_Info, only: Mul
use Definitions, only: iwp

implicit none
integer(kind=iwp), intent(in) :: intentry, isma
integer(kind=iwp), intent(out) :: nlp_value
integer(kind=iwp) :: ijk, ira, ivalue, lra, lri, lritmp, lrj, lrjtmp, lrk, lrktmp, lsmi, lsmij, lsmj, lsmk, m_ia, next_sta, nxo

w0_sdplp25 = (w0_sdplp-w1_sdplp)*w0g25
w1_sdplp25 = (w0_sdplp+w1_sdplp)*w0g25

outer: do LRITMP=NORB_FRZ+1,norb_inn-2
  lsmi = lsm_inn(LRITMP)
  do LRJTMP=LRITMP+1,norb_inn-1
    lsmj = lsm_inn(LRJTMP)
    lsmij = Mul(lsmi,lsmj)
    do LRKTMP=LRJTMP+1,norb_inn
      lsmk = lsm_inn(LRKTMP)
      if (Mul(lsmij,lsmk) /= isma) cycle
      IJK = LRITMP-NORB_FRZ+NGW2(LRJTMP-NORB_FRZ)+NGW3(LRKTMP-NORB_FRZ)
      if (INTIND_IJKA(IJK) == intentry) then
        LRI = LRITMP
        LRJ = LRJTMP
        LRK = LRKTMP
        exit outer
      end if
    end do
  end do
end do outer

next_sta = ibsm_ext(isma)-1

ivalue = 0
do m_ia=1,nlsm_ext(isma)
  ira = m_ia+next_sta
  LRA = norb_number(ira)
  ivalue = ivalue+1
  call TRANS_IJKL_INTPOS(LRA,LRI,LRK,LRJ,NXO)
  index_lpext(ivalue) = NXO
  value_lpext(ivalue) = w0_sdplp25
  call TRANS_IJKL_INTPOS(LRA,LRJ,LRK,LRI,NXO)
  index_lpext1(ivalue) = NXO
  value_lpext1(ivalue) = w1_sdplp25
end do
nlp_value = ivalue

end subroutine lp10_arbrbr_ext_calcuvalue_G

subroutine lp11_arblbr_ext_calcuvalue_G(intentry,isma,nlp_value)

use gugaci_global, only: ibsm_ext, index_lpext, index_lpext1, intind_ijka, lsm_inn, ngw2, ngw3, nlsm_ext, norb_frz, norb_inn, &
                         norb_number, value_lpext, value_lpext1, w0_sdplp, w0_sdplp25, w0g25, w1_sdplp, w1_sdplp25
use Symmetry_Info, only: Mul
use Constants, only: Two
use Definitions, only: iwp

implicit none
integer(kind=iwp), intent(in) :: intentry, isma
integer(kind=iwp), intent(out) :: nlp_value
integer(kind=iwp) :: ijk, ira, ivalue, lra, lri, lritmp, lrj, lrjtmp, lrk, lrktmp, lsmi, lsmij, lsmj, lsmk, m_ia, next_sta, nxo

w0_sdplp25 = (w0_sdplp-w1_sdplp)*w0g25
w1_sdplp25 = -Two*w0_sdplp*w0g25

outer: do LRITMP=NORB_FRZ+1,norb_inn-2
  lsmi = lsm_inn(LRITMP)
  do LRJTMP=LRITMP+1,norb_inn-1
    lsmj = lsm_inn(LRJTMP)
    lsmij = Mul(lsmi,lsmj)
    do LRKTMP=LRJTMP+1,norb_inn
      lsmk = lsm_inn(LRKTMP)
      if (Mul(lsmij,lsmk) /= isma) cycle
      IJK = LRITMP-NORB_FRZ+NGW2(LRJTMP-NORB_FRZ)+NGW3(LRKTMP-NORB_FRZ)
      if (INTIND_IJKA(IJK) == intentry) then
        LRI = LRITMP
        LRJ = LRJTMP
        LRK = LRKTMP
        exit outer
      end if
    end do
  end do
end do outer
next_sta = ibsm_ext(isma)-1

ivalue = 0
do m_ia=1,nlsm_ext(isma)
  ira = m_ia+next_sta
  LRA = norb_number(ira)
  ivalue = ivalue+1
  call TRANS_IJKL_INTPOS(LRA,LRJ,LRK,LRI,NXO)
  index_lpext(ivalue) = NXO
  value_lpext(ivalue) = w0_sdplp25
  call TRANS_IJKL_INTPOS(LRA,LRK,LRJ,LRI,NXO)
  index_lpext1(ivalue) = NXO
  value_lpext1(ivalue) = w1_sdplp25
end do
nlp_value = ivalue

end subroutine lp11_arblbr_ext_calcuvalue_G

subroutine lp12_arblbl_ext_calcuvalue_G(intentry,isma,nlp_value)

use gugaci_global, only: ibsm_ext, index_lpext, index_lpext1, intind_ijka, lsm_inn, ngw2, ngw3, nlsm_ext, norb_frz, norb_inn, &
                         norb_number, value_lpext, value_lpext1, w0_sdplp, w0_sdplp25, w0g25, w1_sdplp, w1_sdplp25
use Symmetry_Info, only: Mul
use Constants, only: Two
use Definitions, only: iwp

implicit none
integer(kind=iwp), intent(in) :: intentry, isma
integer(kind=iwp), intent(out) :: nlp_value
integer(kind=iwp) :: ijk, ira, ivalue, lra, lri, lritmp, lrj, lrjtmp, lrk, lrktmp, lsmi, lsmij, lsmj, lsmk, m_ia, next_sta, nxo

w0_sdplp25 = (w0_sdplp-w1_sdplp)*w0g25
w1_sdplp25 = -Two*w0_sdplp*w0g25

outer: do LRITMP=NORB_FRZ+1,norb_inn-2
  lsmi = lsm_inn(LRITMP)
  do LRJTMP=LRITMP+1,norb_inn-1
    lsmj = lsm_inn(LRJTMP)
    lsmij = Mul(lsmi,lsmj)
    do LRKTMP=LRJTMP+1,norb_inn
      lsmk = lsm_inn(LRKTMP)
      if (Mul(lsmij,lsmk) /= isma) cycle
      IJK = LRITMP-NORB_FRZ+NGW2(LRJTMP-NORB_FRZ)+NGW3(LRKTMP-NORB_FRZ)
      if (INTIND_IJKA(IJK) == intentry) then
        LRI = LRITMP
        LRJ = LRJTMP
        LRK = LRKTMP
        exit outer
      end if
    end do
  end do
end do outer

next_sta = ibsm_ext(isma)-1

ivalue = 0
do m_ia=1,nlsm_ext(isma)
  ira = m_ia+next_sta
  LRA = norb_number(ira)
  ivalue = ivalue+1
  call TRANS_IJKL_INTPOS(LRA,LRI,LRK,LRJ,NXO)
  index_lpext(ivalue) = NXO
  value_lpext(ivalue) = w0_sdplp25
  call TRANS_IJKL_INTPOS(LRA,LRK,LRJ,LRI,NXO)
  index_lpext1(ivalue) = NXO
  value_lpext1(ivalue) = w1_sdplp25
end do
nlp_value = ivalue

end subroutine lp12_arblbl_ext_calcuvalue_G

subroutine lp9_drlbl_ext_calcuvalue_G(lri,lrk,isma)

use gugaci_global, only: ibsm_ext, index_lpext, index_lpext1, nlsm_ext, norb_number, value_lpext, value_lpext1, w0_sdplp, &
                         w0_sdplp25, w0g25, w1_sdplp, w1_sdplp25
use Constants, only: Two
use Definitions, only: iwp

implicit none
integer(kind=iwp), intent(in) :: lri, lrk, isma
integer(kind=iwp) :: ira, ivalue, lra, m_ia, next_sta, nxo

next_sta = ibsm_ext(isma)-1
w0_sdplp25 = (w0_sdplp-w1_sdplp)*w0g25
w1_sdplp25 = -Two*w0_sdplp*w0g25

ivalue = 0

do m_ia=1,nlsm_ext(isma)
  ira = m_ia+next_sta
  LRA = norb_number(ira)
  ivalue = ivalue+1
  call TRANS_IJKL_INTPOS(LRA,LRK,LRI,LRK,NXO)
  index_lpext(ivalue) = NXO
  value_lpext(ivalue) = w0_sdplp25
  call TRANS_IJKL_INTPOS(LRA,LRI,LRK,LRK,NXO)
  index_lpext1(ivalue) = NXO
  value_lpext1(ivalue) = w1_sdplp25
end do

end subroutine lp9_drlbl_ext_calcuvalue_G

subroutine lp8_drlbr_sum_calcuvalue_G(lri,lrk,isma,nv)

use gugaci_global, only: ibsm_ext, index_lpext, index_lpext1, nlsm_ext, norb_number, value_lpext, value_lpext1, w0_sdplp, &
                         w0_sdplp25, w0g25, w1_sdplp25
use Definitions, only: iwp

implicit none
integer(kind=iwp), intent(in) :: lri, lrk, isma
integer(kind=iwp), intent(out) :: nv
integer(kind=iwp) :: ira, ivalue, lra, m_ia, next_sta, nxo

next_sta = ibsm_ext(isma)-1

w0_sdplp25 = w0_sdplp*w0g25
w1_sdplp25 = 2*w0_sdplp*w0g25

ivalue = 0
do m_ia=1,nlsm_ext(isma)
  ira = m_ia+next_sta
  LRA = norb_number(ira)
  ivalue = ivalue+1
  call TRANS_IJKL_INTPOS(LRA,LRK,LRI,LRK,NXO)
  index_lpext(ivalue) = NXO
  value_lpext(ivalue) = w0_sdplp25
  call TRANS_IJKL_INTPOS(LRA,LRI,LRK,LRK,NXO)
  index_lpext1(ivalue) = NXO
  value_lpext1(ivalue) = -w1_sdplp25
end do
nv = ivalue

end subroutine lp8_drlbr_sum_calcuvalue_G

subroutine lp9_drlbl_sum_calcuvalue_G(lri,lrk,isma,nv)

use gugaci_global, only: ibsm_ext, index_lpext, index_lpext1, nlsm_ext, norb_number, value_lpext, value_lpext1, w0_sdplp, &
                         w0_sdplp25, w0g25, w1_sdplp25
use Constants, only: Two
use Definitions, only: iwp

implicit none
integer(kind=iwp), intent(in) :: lri, lrk, isma
integer(kind=iwp), intent(out) :: nv
integer(kind=iwp) :: ira, ivalue, lra, m_ia, next_sta, nxo

next_sta = ibsm_ext(isma)-1
w0_sdplp25 = w0_sdplp*w0g25
w1_sdplp25 = -Two*w0_sdplp*w0g25

ivalue = 0
do m_ia=1,nlsm_ext(isma)
  ira = m_ia+next_sta
  LRA = norb_number(ira)
  ivalue = ivalue+1
  call TRANS_IJKL_INTPOS(LRA,LRK,LRI,LRK,NXO)
  index_lpext(ivalue) = NXO
  value_lpext(ivalue) = w0_sdplp25
  call TRANS_IJKL_INTPOS(LRA,LRI,LRK,LRK,NXO)
  index_lpext1(ivalue) = NXO
  value_lpext1(ivalue) = w1_sdplp25
end do
nv = ivalue

end subroutine lp9_drlbl_sum_calcuvalue_G

subroutine gsd_ext_sequence_G(iltype,ilsm,irsm,lri)

use gugaci_global, only: ibsm_ext, icano_nnend, icano_nnsta, icnt_base, iesm_ext, iseg_downwei, isegdownwei, m_jc, m_jd, ng_sm
use Symmetry_Info, only: Mul
use Definitions, only: iwp

implicit none
integer(kind=iwp), intent(in) :: iltype, ilsm, irsm, lri
integer(kind=iwp) :: ic, icano_nn, icend, icsta, ilnodedownwei, indl, isma, ismb, ismnoded, ismnodes

ismnodes = ilsm
ismnoded = irsm
indl = 0 !?
if (iltype == 2) indl = 1+ismnodes
if (iltype == 3) indl = 9+ismnodes
if (iltype == 4) indl = 17+ismnodes
ilnodedownwei = iseg_downwei(indl)
isegdownwei = ilnodedownwei
icano_nnsta = 1
icnt_base = 0
icsta = ibsm_ext(ismnoded)
icend = iesm_ext(ismnoded)
m_jc = 0

do ic=icsta,icend
  m_jd = ic
  m_jc = ic-icsta+1
  icano_nn = m_jc
  icano_nnend = icano_nn
  do ismb=1,ismnoded-1
    isma = Mul(ismnodes,ismb)
    if (isma > ismb) cycle
    call g31_diffsym_G(lri,isma,ismb)
  end do

  ismb = ismnoded
  isma = Mul(ismnodes,ismb)
  if (isma == ismb) then
    call gsd_samesym_aaa_G(lri,isma)
  else if (isma < ismb) then
    call gsd_diffsamesym_abb_G(lri,isma,ismb)
  end if

  do ismb=ismnoded+1,ng_sm
    isma = Mul(ismnodes,ismb)
    if (isma > ismb) cycle
    if (ismnoded > isma) then
      call g32a_diffsym_G(lri,isma,ismb)
    else if (ismnoded == isma) then
      call gsd_diffsamesym_aab_G(lri,isma,ismb)
    else
      call g32b_diffsym_G(lri,isma,ismb)
    end if
  end do

  if ((ismnodes == 1) .and. (iltype == 4)) then
    call gsd_arlp_s1_G(lri)
  end if
  icnt_base = icnt_base+ilnodedownwei
end do

end subroutine gsd_ext_sequence_G

subroutine lp678_ext_calcuvalue_G(lri,lrk,isma,nlp_value)

use gugaci_global, only: ibsm_ext, iesm_ext, index_lpext, index_lpext1, nlsm_ext, norb_number, value_lpext, w0_sdplp, w0_sdplp25, &
                         w0g25
use Definitions, only: iwp

implicit none
integer(kind=iwp), intent(in) :: lri, lrk, isma
integer(kind=iwp), intent(out) :: nlp_value
integer(kind=iwp) :: ia, iaend, iasta, ilpvalue, lra, nxo

w0_sdplp25 = w0_sdplp*w0g25
ilpvalue = 0
iasta = ibsm_ext(isma)
iaend = iesm_ext(isma)
do ia=iasta,iaend
  LRA = norb_number(ia)
  ilpvalue = ilpvalue+1
  call TRANS_IJKL_INTPOS(LRA,LRK,LRI,LRK,NXO)
  index_lpext(ilpvalue) = NXO
  value_lpext(ilpvalue) = w0_sdplp25
  index_lpext1(ilpvalue) = 0
end do
nlp_value = nlsm_ext(isma)

end subroutine lp678_ext_calcuvalue_G

subroutine lp_ar_coe_calcuvalue_G(idtu,isma,lri,lrj,nlp_value,lpcoe,nlp_value1)

use gugaci_global, only: ibsm_ext, ican_a, index_lpext3, index_lpext4, index_lpext5, jb_sys, logic_dh, nlsm_ext, norb_dz, &
                         norb_inn, norb_number, value_lpext3, value_lpext4, value_lpext5, w0_sdplp, w0_sdplp25, w0g25
use Constants, only: Two
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: idtu, isma, lri, lrj, lpcoe(norb_dz+1:norb_inn)
integer(kind=iwp), intent(out) :: nlp_value, nlp_value1
integer(kind=iwp) :: icoe, idorb, ilpvalue, ilpvalue1, iorb, iorbs, ira, kcoe, lend, lra, lsta, m_ia, ndorb, next_sta, nia, nocc, &
                     nsorb, nxo
real(kind=wp) :: tcoe

next_sta = ibsm_ext(isma)-1
w0_sdplp25 = w0_sdplp*w0g25
ilpvalue = 0

if (idtu == 100) then
  lsta = lri
  lend = norb_inn
  ilpvalue1 = 0
  do m_ia=1,nlsm_ext(isma)
    ira = m_ia+next_sta
    lra = norb_number(ira)
    NIA = ICAN_A(LRA)+LRI
    ilpvalue = ilpvalue+1
    index_lpext5(ilpvalue) = NIA
    value_lpext5(ilpvalue) = w0_sdplp25

    ilpvalue1 = 0
    do iorb=lsta,lend
      KCOE = lpcoe(iorb)
      call NEOC(KCOE,NOCC,TCOE)

      ilpvalue1 = ilpvalue1+1
      call TRANS_IJKL_INTPOS(LRA,iorb,LRI,iorb,NXO)
      index_lpext3(ilpvalue,ilpvalue1) = NXO
      value_lpext3(ilpvalue,ilpvalue1) = w0_sdplp25*nocc*tcoe
      call TRANS_IJKL_INTPOS(LRA,LRI,iorb,iorb,NXO)
      index_lpext4(ilpvalue,ilpvalue1) = NXO
      value_lpext4(ilpvalue,ilpvalue1) = w0_sdplp25*nocc

    end do
  end do
  nlp_value = nlsm_ext(isma)
  nlp_value1 = ilpvalue1
  return
end if

if (idtu == 51) then
  lsta = lri
  lend = norb_DZ
  ilpvalue1 = 0
  do m_ia=1,nlsm_ext(isma)
    ira = m_ia+next_sta
    lra = norb_number(ira)
    NIA = ICAN_A(LRA)+LRI
    ilpvalue = ilpvalue+1
    index_lpext5(ilpvalue) = NIA
    value_lpext5(ilpvalue) = w0_sdplp25

    ilpvalue1 = 0
    if (LOGIC_DH) then
      ilpvalue1 = ilpvalue1+1
      call TRANS_IJKL_INTPOS(LRA,LRI,LRI,LRI,NXO)
      index_lpext3(ilpvalue,ilpvalue1) = NXO
      value_lpext3(ilpvalue,ilpvalue1) = w0_sdplp25
      index_lpext4(ilpvalue,ilpvalue1) = 0

      do idorb=lsta+1,lend
        ilpvalue1 = ilpvalue1+1
        call TRANS_IJKL_INTPOS(LRA,idorb,LRI,idorb,NXO)
        index_lpext3(ilpvalue,ilpvalue1) = NXO
        value_lpext3(ilpvalue,ilpvalue1) = -w0_sdplp25
        call TRANS_IJKL_INTPOS(LRA,LRI,idorb,idorb,NXO)
        index_lpext4(ilpvalue,ilpvalue1) = NXO
        value_lpext4(ilpvalue,ilpvalue1) = w0_sdplp25*Two

      end do
    end if
    iorbs = max(norb_dz+1,lri)
    do iorb=iorbs,norb_inn
      Kcoe = lpcoe(iorb)
      call NEOC(KCOE,NOCC,TCOE)
      ilpvalue1 = ilpvalue1+1
      call TRANS_IJKL_INTPOS(LRA,iorb,LRI,iorb,NXO)
      index_lpext3(ilpvalue,ilpvalue1) = NXO
      value_lpext3(ilpvalue,ilpvalue1) = w0_sdplp25*nocc*tcoe
      call TRANS_IJKL_INTPOS(LRA,LRI,iorb,iorb,NXO)
      index_lpext4(ilpvalue,ilpvalue1) = NXO
      value_lpext4(ilpvalue,ilpvalue1) = w0_sdplp25*nocc
    end do
  end do
  nlp_value = nlsm_ext(isma)
  nlp_value1 = ilpvalue1
  return
end if

if ((idtu == 25) .or. (idtu == 43)) then
  ilpvalue1 = 0
  do m_ia=1,nlsm_ext(isma)
    ira = m_ia+next_sta
    lra = norb_number(ira)
    NIA = ICAN_A(LRA)+LRI
    ilpvalue = ilpvalue+1
    index_lpext5(ilpvalue) = NIA
    value_lpext5(ilpvalue) = w0_sdplp25

    ilpvalue1 = 0

    ilpvalue1 = ilpvalue1+1
    call TRANS_IJKL_INTPOS(LRA,LRI,LRI,LRI,NXO)
    index_lpext3(ilpvalue,ilpvalue1) = NXO
    value_lpext3(ilpvalue,ilpvalue1) = w0_sdplp25
    index_lpext4(ilpvalue,ilpvalue1) = 0

    do idorb=lri+1,norb_dz
      ilpvalue1 = ilpvalue1+1
      call TRANS_IJKL_INTPOS(LRA,idorb,LRI,idorb,NXO)
      index_lpext3(ilpvalue,ilpvalue1) = NXO
      value_lpext3(ilpvalue,ilpvalue1) = -w0_sdplp25
      call TRANS_IJKL_INTPOS(LRA,LRI,idorb,idorb,NXO)
      index_lpext4(ilpvalue,ilpvalue1) = NXO
      value_lpext4(ilpvalue,ilpvalue1) = w0_sdplp25*Two
    end do
    iorbs = max(norb_dz+1,lri)
    do iorb=iorbs,norb_inn
      Kcoe = lpcoe(iorb)
      call NEOC(KCOE,NOCC,TCOE)
      ilpvalue1 = ilpvalue1+1
      call TRANS_IJKL_INTPOS(LRA,iorb,LRI,iorb,NXO)
      index_lpext3(ilpvalue,ilpvalue1) = NXO
      value_lpext3(ilpvalue,ilpvalue1) = w0_sdplp25*nocc*tcoe
      call TRANS_IJKL_INTPOS(LRA,LRI,iorb,iorb,NXO)
      index_lpext4(ilpvalue,ilpvalue1) = NXO
      value_lpext4(ilpvalue,ilpvalue1) = w0_sdplp25*nocc
    end do
  end do
  nlp_value = nlsm_ext(isma)
  nlp_value1 = ilpvalue1
  return
end if

if (idtu == 26) then
  lsta = lri
  lend = norb_DZ
  ilpvalue1 = 0
  do m_ia=1,nlsm_ext(isma)
    ira = m_ia+next_sta
    lra = norb_number(ira)
    NIA = ICAN_A(LRA)+LRI
    ilpvalue = ilpvalue+1
    index_lpext5(ilpvalue) = NIA
    value_lpext5(ilpvalue) = w0_sdplp25

    ilpvalue1 = 0

    if (LOGIC_DH) then
      do idorb=lsta+1,lend
        ilpvalue1 = ilpvalue1+1
        call TRANS_IJKL_INTPOS(LRA,idorb,LRI,idorb,NXO)
        index_lpext3(ilpvalue,ilpvalue1) = NXO
        value_lpext3(ilpvalue,ilpvalue1) = -w0_sdplp25
        call TRANS_IJKL_INTPOS(LRA,LRI,idorb,idorb,NXO)
        index_lpext4(ilpvalue,ilpvalue1) = NXO
        value_lpext4(ilpvalue,ilpvalue1) = w0_sdplp25*Two
      end do
    end if
    do iorb=norb_dz+1,norb_inn
      Kcoe = lpcoe(iorb)
      call NEOC(KCOE,NOCC,TCOE)
      ilpvalue1 = ilpvalue1+1
      call TRANS_IJKL_INTPOS(LRA,iorb,LRI,iorb,NXO)
      index_lpext3(ilpvalue,ilpvalue1) = NXO
      value_lpext3(ilpvalue,ilpvalue1) = w0_sdplp25*nocc*tcoe
      call TRANS_IJKL_INTPOS(LRA,LRI,iorb,iorb,NXO)
      index_lpext4(ilpvalue,ilpvalue1) = NXO
      value_lpext4(ilpvalue,ilpvalue1) = w0_sdplp25*nocc
    end do
  end do
  nlp_value = nlsm_ext(isma)
  nlp_value1 = ilpvalue1
  return
end if
if ((idtu == 28) .or. (idtu == 46)) then
  nsorb = 1
  ndorb = norb_DZ-lri-1
  icoe = 0
  if (ndorb > 0) nsorb = 1
  if (idtu == 28) icoe = -(JB_SYS+2)
  if (idtu == 46) icoe = 0
  ilpvalue1 = 0
  do m_ia=1,nlsm_ext(isma)
    ira = m_ia+next_sta
    lra = norb_number(ira)
    NIA = ICAN_A(LRA)+LRI
    ilpvalue = ilpvalue+1
    index_lpext5(ilpvalue) = NIA
    value_lpext5(ilpvalue) = w0_sdplp25

    ilpvalue1 = 0

    ilpvalue1 = ilpvalue1+1
    call TRANS_IJKL_INTPOS(LRA,LRI,LRI,LRI,NXO)
    index_lpext3(ilpvalue,ilpvalue1) = NXO
    value_lpext3(ilpvalue,ilpvalue1) = w0_sdplp25
    index_lpext4(ilpvalue,ilpvalue1) = 0

    if (nsorb > 0) then
      do iorb=lri+1,norb_dz
        if (iorb == lrj) then
          ilpvalue1 = ilpvalue1+1
          call TRANS_IJKL_INTPOS(LRA,iorb,LRI,iorb,NXO)
          index_lpext3(ilpvalue,ilpvalue1) = NXO
          value_lpext3(ilpvalue,ilpvalue1) = w0_sdplp25*icoe
          call TRANS_IJKL_INTPOS(LRA,LRI,iorb,iorb,NXO)
          index_lpext4(ilpvalue,ilpvalue1) = NXO
          value_lpext4(ilpvalue,ilpvalue1) = w0_sdplp25
        end if
        if (iorb /= lrj) then
          ilpvalue1 = ilpvalue1+1
          call TRANS_IJKL_INTPOS(LRA,iorb,LRI,iorb,NXO)
          index_lpext3(ilpvalue,ilpvalue1) = NXO
          value_lpext3(ilpvalue,ilpvalue1) = -w0_sdplp25
          call TRANS_IJKL_INTPOS(LRA,LRI,iorb,iorb,NXO)
          index_lpext4(ilpvalue,ilpvalue1) = NXO
          value_lpext4(ilpvalue,ilpvalue1) = w0_sdplp25*Two
        end if
      end do
    end if
    do iorb=norb_dz+1,norb_inn
      Kcoe = lpcoe(iorb)
      call NEOC(KCOE,NOCC,TCOE)
      ilpvalue1 = ilpvalue1+1
      call TRANS_IJKL_INTPOS(LRA,iorb,LRI,iorb,NXO)
      index_lpext3(ilpvalue,ilpvalue1) = NXO
      value_lpext3(ilpvalue,ilpvalue1) = w0_sdplp25*nocc*tcoe
      call TRANS_IJKL_INTPOS(LRA,LRI,iorb,iorb,NXO)
      index_lpext4(ilpvalue,ilpvalue1) = NXO
      value_lpext4(ilpvalue,ilpvalue1) = w0_sdplp25*nocc
    end do
  end do
  nlp_value = nlsm_ext(isma)
  nlp_value1 = ilpvalue1
  return
end if
if ((idtu == 57) .or. (idtu == 29)) then
  nsorb = 1
  ndorb = norb_DZ-lri-1
  if (ndorb > 0) nsorb = 1
  if (idtu == 29) then
    icoe = JB_SYS
  else ! if (idtu == 57)
    icoe = -1
  end if
  ilpvalue1 = 0
  do m_ia=1,nlsm_ext(isma)
    ira = m_ia+next_sta
    lra = norb_number(ira)
    NIA = ICAN_A(LRA)+LRI
    ilpvalue = ilpvalue+1
    index_lpext5(ilpvalue) = NIA
    value_lpext5(ilpvalue) = w0_sdplp25

    ilpvalue1 = 0

    ilpvalue1 = ilpvalue1+1
    call TRANS_IJKL_INTPOS(LRA,LRI,LRI,LRI,NXO)
    index_lpext3(ilpvalue,ilpvalue1) = NXO
    value_lpext3(ilpvalue,ilpvalue1) = w0_sdplp25
    index_lpext4(ilpvalue,ilpvalue1) = 0

    if (nsorb > 0) then
      do iorb=lri+1,norb_dz
        if (iorb == lrj) then
          ilpvalue1 = ilpvalue1+1
          call TRANS_IJKL_INTPOS(LRA,iorb,LRI,iorb,NXO)
          index_lpext3(ilpvalue,ilpvalue1) = NXO
          value_lpext3(ilpvalue,ilpvalue1) = w0_sdplp25*icoe
          call TRANS_IJKL_INTPOS(LRA,LRI,iorb,iorb,NXO)
          index_lpext4(ilpvalue,ilpvalue1) = NXO
          value_lpext4(ilpvalue,ilpvalue1) = w0_sdplp25
        end if
        if (iorb /= lrj) then
          ilpvalue1 = ilpvalue1+1
          call TRANS_IJKL_INTPOS(LRA,iorb,LRI,iorb,NXO)
          index_lpext3(ilpvalue,ilpvalue1) = NXO
          value_lpext3(ilpvalue,ilpvalue1) = -w0_sdplp25
          call TRANS_IJKL_INTPOS(LRA,LRI,iorb,iorb,NXO)
          index_lpext4(ilpvalue,ilpvalue1) = NXO
          value_lpext4(ilpvalue,ilpvalue1) = w0_sdplp25*Two
        end if
      end do
    end if

    do iorb=norb_dz+1,norb_inn
      Kcoe = lpcoe(iorb)
      call NEOC(KCOE,NOCC,TCOE)
      ilpvalue1 = ilpvalue1+1
      call TRANS_IJKL_INTPOS(LRA,iorb,LRI,iorb,NXO)
      index_lpext3(ilpvalue,ilpvalue1) = NXO
      value_lpext3(ilpvalue,ilpvalue1) = w0_sdplp25*nocc*tcoe
      call TRANS_IJKL_INTPOS(LRA,LRI,iorb,iorb,NXO)
      index_lpext4(ilpvalue,ilpvalue1) = NXO
      value_lpext4(ilpvalue,ilpvalue1) = w0_sdplp25*nocc
    end do
  end do
  nlp_value = nlsm_ext(isma)
  nlp_value1 = ilpvalue1
  return
end if
!=======================================================================
if ((idtu == 55) .or. (idtu == 73)) then    !(11)(23)=55 (11)(13)=75
  ilpvalue1 = 0
  do m_ia=1,nlsm_ext(isma)
    ira = m_ia+next_sta
    lra = norb_number(ira)
    NIA = ICAN_A(LRA)+LRI
    ilpvalue = ilpvalue+1
    index_lpext5(ilpvalue) = NIA
    value_lpext5(ilpvalue) = w0_sdplp25

    ilpvalue1 = 0

    ilpvalue1 = ilpvalue1+1
    call TRANS_IJKL_INTPOS(LRA,LRI,LRI,LRI,NXO)
    index_lpext3(ilpvalue,ilpvalue1) = NXO
    value_lpext3(ilpvalue,ilpvalue1) = w0_sdplp25
    index_lpext4(ilpvalue,ilpvalue1) = 0
    do idorb=lri+1,norb_dz
      ilpvalue1 = ilpvalue1+1
      call TRANS_IJKL_INTPOS(LRA,idorb,LRI,idorb,NXO)
      index_lpext3(ilpvalue,ilpvalue1) = NXO
      value_lpext3(ilpvalue,ilpvalue1) = -w0_sdplp25
      call TRANS_IJKL_INTPOS(LRA,LRI,idorb,idorb,NXO)
      index_lpext4(ilpvalue,ilpvalue1) = NXO
      value_lpext4(ilpvalue,ilpvalue1) = w0_sdplp25*Two
    end do
    iorbs = max(norb_dz+1,lri)
    do iorb=iorbs,norb_inn
      Kcoe = lpcoe(iorb)
      call NEOC(KCOE,NOCC,TCOE)
      ilpvalue1 = ilpvalue1+1
      call TRANS_IJKL_INTPOS(LRA,iorb,LRI,iorb,NXO)
      index_lpext3(ilpvalue,ilpvalue1) = NXO
      value_lpext3(ilpvalue,ilpvalue1) = w0_sdplp25*nocc*tcoe
      call TRANS_IJKL_INTPOS(LRA,LRI,iorb,iorb,NXO)
      index_lpext4(ilpvalue,ilpvalue1) = NXO
      value_lpext4(ilpvalue,ilpvalue1) = w0_sdplp25*nocc
    end do
  end do
  nlp_value = nlsm_ext(isma)
  nlp_value1 = ilpvalue1
  return
end if

end subroutine lp_ar_coe_calcuvalue_G

subroutine lp_arbl_ext_dd_calcuvalue_G(lri,lrj,iml,imr,nlp_value)

use gugaci_global, only: ibsm_ext, iesm_ext, index_lpext, index_lpext1, int_dd_drl, logic_g49a, logic_g49b, logic_g50, &
                         norb_number, value_lpext, value_lpext1, w0_plp, w0gdd, w1_plp, w1gdd
use Constants, only: Two
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: lri, lrj, iml, imr
integer(kind=iwp), intent(out) :: nlp_value
integer(kind=iwp) :: ira, irb, ivalue, lra, lrb, nlbf, nlef, nrbf, nref, nxo
real(kind=wp) :: valuetmp1, w0lp, w1lp

nlbf = ibsm_ext(iml)
nlef = iesm_ext(iml)
nrbf = ibsm_ext(imr)
nref = iesm_ext(imr)

w0lp = w0_plp*w0gdd
w1lp = w1_plp*w1gdd
valuetmp1 = w0lp
w0lp = w0lp-w1lp
w1lp = -valuetmp1*Two
ivalue = 0
! G50
if (logic_g50) then
  if (logic_g49b) then

    do ira=nlbf,nlef
      lra = norb_number(ira)
      ivalue = ivalue+1
      call TRANS_IJKL_INTPOS(lrj,lra,lri,lra,NXO)
      index_lpext(ivalue) = NXO
      value_lpext(ivalue) = w0lp
      call TRANS_IJKL_INTPOS(lrj,lri,lra,lra,NXO)
      index_lpext1(ivalue) = NXO
      value_lpext1(ivalue) = w1lp
    end do
  end if

  ivalue = ivalue+int_dd_drl

  do irb=nlbf,nlef
    lrb = norb_number(irb)
    !lsmb = lsm(irb)
    do ira=nlbf,irb-1
      lra = norb_number(ira)
      !lsma = lsm(ira)
      !lsmba = Mul(lsmb,lsma)
      ivalue = ivalue+1
      call TRANS_IJKL_INTPOS(lra,lri,lrj,lrb,NXO)
      index_lpext(ivalue) = NXO
      value_lpext(ivalue) = w0lp
      call TRANS_IJKL_INTPOS(lra,lrb,lrj,lri,NXO)
      index_lpext1(ivalue) = NXO
      value_lpext1(ivalue) = w1lp
    end do
  end do
  do irb=nlbf,nlef
    lrb = norb_number(irb)
    !lsmb = lsm(irb)
    do ira=nlbf,irb-1
      lra = norb_number(ira)
      !lsma = lsm(ira)
      !lsmba = Mul(lsmb,lsma)
      ivalue = ivalue+1
      call TRANS_IJKL_INTPOS(lra,lrj,lrb,lri,NXO)
      index_lpext(ivalue) = NXO
      value_lpext(ivalue) = w0lp
      call TRANS_IJKL_INTPOS(lra,lrb,lrj,lri,NXO)
      index_lpext1(ivalue) = NXO
      value_lpext1(ivalue) = w1lp
    end do
  end do

else
  ivalue = ivalue+int_dd_drl
  if (logic_g49a) then
    ! G49a:Bl_Ar  line=12
    do irb=nrbf,nref
      lrb = norb_number(irb)
      do ira=nlbf,nlef
        lra = norb_number(ira)

        ivalue = ivalue+1
        call TRANS_IJKL_INTPOS(lra,lri,lrj,lrb,NXO)
        index_lpext(ivalue) = NXO
        value_lpext(ivalue) = w0lp
        call TRANS_IJKL_INTPOS(lra,lrb,lrj,lri,NXO)
        index_lpext1(ivalue) = NXO
        value_lpext1(ivalue) = w1lp
      end do
    end do
  else
    ! G49b:Br_Al  line=11
    do irb=nlbf,nlef
      lrb = norb_number(irb)
      do ira=nrbf,nref
        lra = norb_number(ira)

        ivalue = ivalue+1
        call TRANS_IJKL_INTPOS(lra,lrj,lrb,lri,NXO)
        index_lpext(ivalue) = NXO
        value_lpext(ivalue) = w0lp
        call TRANS_IJKL_INTPOS(lra,lrb,lrj,lri,NXO)
        index_lpext1(ivalue) = NXO
        value_lpext1(ivalue) = w1lp
      end do
    end do
  end if
end if
nlp_value = ivalue

end subroutine lp_arbl_ext_dd_calcuvalue_G

subroutine lp_drl_ext_dd_calcuvalue_G(lri,iml,nlp_value)

use gugaci_global, only: ibsm_ext, iesm_ext, index_lpext, index_lpext1, int_dd_drl, logic_g49b, nlsm_ext, norb_number, &
                         value_lpext, value_lpext1, w0_plp, w0gdd, w1_plp, w1gdd
use Constants, only: Two
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: lri, iml
integer(kind=iwp), intent(out) :: nlp_value
integer(kind=iwp) :: iaend, iasta, ira, irb, ivalue, jvalue, lra, lrb, mloop, nliml, nxo
real(kind=wp) :: w0lp, w1lp

w0lp = w0_plp*w0gdd
w1lp = w1_plp*w1gdd
nliml = nlsm_ext(iml)
iasta = ibsm_ext(iml)
iaend = iesm_ext(iml)

ivalue = 0
if (logic_g49b) then
  do ira=iasta,iaend
    lra = norb_number(ira)
    ivalue = ivalue+1
    call TRANS_IJKL_INTPOS(lra,LRI,lra,LRI,NXO)
    index_lpext(ivalue) = NXO
    value_lpext(ivalue) = -w1lp*Two
    index_lpext1(ivalue) = 0
  end do
end if
mloop = nliml*(nliml-1)/2

ivalue = ivalue+int_dd_drl
jvalue = 0
do irb=iasta,iaend
  lrb = norb_number(irb)
  do ira=iasta,irb-1
    lra = norb_number(ira)
    ivalue = ivalue+1
    jvalue = ivalue+mloop
    call TRANS_IJKL_INTPOS(lra,LRI,lrb,LRI,NXO)
    index_lpext(ivalue) = NXO
    value_lpext(ivalue) = w0lp-w1lp
    index_lpext(jvalue) = index_lpext(ivalue)
    value_lpext(jvalue) = value_lpext(ivalue)

    call TRANS_IJKL_INTPOS(lra,lrb,LRI,LRI,NXO)
    index_lpext1(ivalue) = NXO
    value_lpext1(ivalue) = -Two*w0lp
    index_lpext1(jvalue) = index_lpext1(ivalue)
    value_lpext1(jvalue) = value_lpext1(ivalue)
  end do
end do
nlp_value = jvalue

end subroutine lp_drl_ext_dd_calcuvalue_G

subroutine lp_arbr_ext_svtv_calcuvalue_G(lri,lrj,nlp_value)

use gugaci_global, only: ibsm_ext, iesm_ext, index_lpext, index_lpext1, logic_g13, lsm_inn, ng_sm, norb_ext, norb_number, &
                         value_lpext, value_lpext1, w0_plp, w0g13a, w0g36a, w1_plp, w1g36a
use Symmetry_Info, only: Mul
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: lri, lrj
integer(kind=iwp), intent(out) :: nlp_value
integer(kind=iwp) :: ic, icend, icsta, id, idend, idsta, ivalue, jc, jd, lrc, lsmc, lsmd, lsmi, lsmij, lsmj, nxo
real(kind=wp) :: valuelptmp1, w0lp, w1lp

ivalue = 0
lsmi = lsm_inn(lri)
lsmj = lsm_inn(lrj)
lsmij = Mul(lsmi,lsmj)
! G36a
w0lp = w0_plp*w0g36a
w1lp = w1_plp*w1g36a
valuelptmp1 = w0lp
w0lp = w0lp-w1lp
w1lp = valuelptmp1+w1lp
! ArBr -- B^rA^r =10
do lsmc=1,ng_sm
  lsmd = Mul(lsmij,lsmc)
  if (lsmd > lsmc) cycle
  icsta = ibsm_ext(lsmc)
  icend = iesm_ext(lsmc)
  idsta = ibsm_ext(lsmd)
  idend = iesm_ext(lsmd)
  if (lsmc == lsmd) icsta = icsta+1
  do ic=icsta,icend
    jc = norb_number(ic)
    do id=idsta,min(idend,ic-1)
      jd = norb_number(id)
      ivalue = ivalue+1
      call TRANS_IJKL_INTPOS(jd,lri,lrj,jc,NXO)
      index_lpext(ivalue) = NXO
      value_lpext(ivalue) = w0lp
      call TRANS_IJKL_INTPOS(jd,lrj,jc,lri,NXO)
      index_lpext1(ivalue) = NXO
      value_lpext1(ivalue) = w1lp
    end do
  end do
end do
! G36b
! G1415
if (logic_g13) then

  w0lp = w0g13a*(w0_plp+w1_plp)
  do ic=1,norb_ext
    lrc = norb_number(ic)
    ivalue = ivalue+1
    call TRANS_IJKL_INTPOS(lrj,lrc,lri,lrc,NXO)
    index_lpext(ivalue) = NXO
    value_lpext(ivalue) = w0lp
    index_lpext1(ivalue) = 0
  end do
end if
nlp_value = ivalue

end subroutine lp_arbr_ext_svtv_calcuvalue_G

subroutine lp_drr_ext_svtv_calcuvalue_G(lri,nlp_value)

use gugaci_global, only: ibsm_ext, iesm_ext, index_lpext, index_lpext1, logic_g13, ng_sm, norb_all, norb_inn, norb_number, &
                         value_lpext, w0_plp, w0g13a, w0g36a, w1_plp
use Constants, only: Two, Half
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: lri
integer(kind=iwp), intent(out) :: nlp_value
integer(kind=iwp) :: ibend, ibsta, ira, irb, ivalue, lmb, lra, lrb, nxo
real(kind=wp) :: w0lp

ivalue = 0
! G36a
w0lp = (w0_plp+w1_plp)*w0g36a
do lmb=1,ng_sm
  ibsta = ibsm_ext(lmb)
  ibend = iesm_ext(lmb)
  do irb=ibsta,ibend
    lrb = norb_number(irb)
    do ira=ibsta,irb-1
      lra = norb_number(ira)
      ivalue = ivalue+1
      call TRANS_IJKL_INTPOS(lra,lri,lrb,lri,NXO)
      index_lpext(ivalue) = NXO
      value_lpext(ivalue) = w0lp
      index_lpext1(ivalue) = 0
    end do
  end do
end do

if (logic_g13) then
!=======================================================================
! Drr-DRR
! w0lp=Two*w0lp but not One*w0lp is based on that the non-diagonal
! just uses the non-triangle <Ci|H|Cj> which designates that I > J.

  w0lp = Half*w0_plp*w0g13a
  w0lp = Two*w0lp
  do LRA=NORB_ALL,NORB_INN+1,-1
    ivalue = ivalue+1
    call TRANS_IJKL_INTPOS(LRA,lri,LRA,lri,NXO)
    index_lpext(ivalue) = NXO
    value_lpext(ivalue) = w0lp
    index_lpext1(ivalue) = 0
  end do
end if
nlp_value = ivalue

end subroutine lp_drr_ext_svtv_calcuvalue_G
