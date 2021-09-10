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

! SUBROUTINE lp10_arbrbr_ext_calcuvalue(intentry,isma,nlp_value)
! SUBROUTINE lp11_arblbr_ext_calcuvalue(intentry,isma,nlp_value)
! SUBROUTINE lp12_arblbl_ext_calcuvalue(intentry,isma,nlp_value)
! SUBROUTINE lp9_drlbl_ext_calcuvalue(lri,lrk,isma)
! SUBROUTINE lp8_drlbr_sum_calcuvalue(lri,lrp,lrq,isma,nv)
! SUBROUTINE lp9_drlbl_sum_calcuvalue(lri,lrp,lrq,isma,nv)
! SUBROUTINE lp_arbl_ext_dd_calcuvalue(lri,lrj,iml,imr,nlp_value)
! SUBROUTINE lp_ar_coe_calcuvalue(idtu,isma,lri,lrj,nlp_value,lpcoe
! SUBROUTINE lp_drl_ext_SS_calcuvalue(lri,nlp_value)
! SUBROUTINE lp_drl_sum_SS_calcuvalue(lri,lrj,nlp_value)
! SUBROUTINE lp_drl_ext_ST_calcuvalue(lri,nlp_value)
! SUBROUTINE lp_drl_ext_TT_calcuvalue(lri,n1415_value,nlp_value)
! SUBROUTINE lp_drl_sum_TT_calcuvalue(lri,lrj,n1415,nlp_value)
! SUBROUTINE lp_arbr_ext_svtv_calcuvalue(LRI,LRJ,nlp_value)
! SUBROUTINE lp_drr_ext_svtv_calcuvalue(lri,nlp_value)
! SUBROUTINE lp_arbl_ext_st_calcuvalue(lri,lrj,nlp_value)
! SUBROUTINE lp678_ext_calcuvalue(lri,lrk,isma,nlp_value)

subroutine lp_drl_ext_TS_calcuvalue(lri,nlp_value)

#include "drt_h.fh"
#include "lpextmode_h.fh"
#include "intsort_h.fh"

intpos = intind_abkk(lri)
intspace = intspace_abkk(lri)
ivalue = 0
! G1415
if (logic_g1415) then

  w1lp = w1_plp*w1g14a

  do ismb=1,ng_sm
    isma = mul_tab(ismb,ism_g1415)
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
        !value_lpext(ivalue) = vint_ci(intposia)*ww0lp+vint_ci(intposia+1)*ww1lp+valuelpib
        ! OK,only for Spin=0
        value_lpext(ivalue) = (voint(lra,lri)-voint(lrb,lri))*w1lp
      end do
    end do
  end do
end if

! G2G4b
if (logic_g2g4a) then
  w0lp = w0_plp*w0g2a
  w1lp = w1_plp*w1g2a

! valuelptmp1 = w0lp
! w0lp = w0lp-w1lp
! w1lp = -valuelptmp1*2.d0
! valuelptmp1 = ww0lp
! ww0lp = ww0lp-ww1lp
! ww1lp = -valuelptmp1*2.d0

  do i=1,intspace
    ivalue = ivalue+2
    ! Drl -- B^rA^l =4_2
    value_lpext(ivalue) = w1lp*vint_ci(intpos)
    ! Drl -- B^lA^r =4_3
    value_lpext(ivalue-1) = -value_lpext(ivalue)
    intpos = intpos+2
  end do
end if

intpos = intind_abkk(lri)
w0lp = w0_plp*w0g36a
w1lp = w1_plp*w1g36a

do I=1,intspace
  ivalue = ivalue+1
  value_lpext(ivalue) = vint_ci(intpos+1)*w0lp-vint_ci(intpos)*w1lp
  intpos = intpos+2
end do
nlp_value = ivalue

end subroutine lp_drl_ext_TS_calcuvalue

subroutine lp9_drlbl_ext_sd_calcuvalue(intentry,isma)   ! ,nlp_value)

#include "drt_h.fh"
#include "intsort_h.fh"
#include "lpextmode_h.fh"

w0_sdplp25 = (w0_sdplp-w1_sdplp)*w0g25
w1_sdplp25 = -2.d0*w0_sdplp*w0g25
intpos = intentry
iaddpos = norb_act*2            !severe_new_error_1020
ilpvalue = 0
do m_ia=1,nlsm_ext(isma)
  ilpvalue = ilpvalue+1
  value_lpext(ilpvalue) = w0_sdplp25*vint_ci(intpos)+w1_sdplp25*vint_ci(intpos+1)
  intpos = intpos+iaddpos      !ip3ad_intspace   !severe_new_error_
end do

end subroutine lp9_drlbl_ext_sd_calcuvalue

subroutine lp10_arbrbr_ext_calcuvalue(intentry,isma,nlp_value)

#include "drt_h.fh"
#include "intsort_h.fh"
#include "lpextmode_h.fh"

w0_sdplp25 = (w0_sdplp-w1_sdplp)*w0g25
w1_sdplp25 = (w0_sdplp+w1_sdplp)*w0g25
intpos = intentry
ilpvalue = 0
do m_ia=1,nlsm_ext(isma)
  ilpvalue = ilpvalue+1
  value_lpext(ilpvalue) = w0_sdplp25*vint_ci(intpos+2)+w1_sdplp25*vint_ci(intpos+1)
  intpos = intpos+3
end do
nlp_value = ilpvalue

end subroutine lp10_arbrbr_ext_calcuvalue

subroutine lp11_arblbr_ext_calcuvalue(intentry,isma,nlp_value)

#include "drt_h.fh"
#include "intsort_h.fh"
#include "lpextmode_h.fh"

w0_sdplp25 = (w0_sdplp-w1_sdplp)*w0g25
w1_sdplp25 = -2.d0*w0_sdplp*w0g25
intpos = intentry
ilpvalue = 0
do m_ia=1,nlsm_ext(isma)
  ilpvalue = ilpvalue+1
  value_lpext(ilpvalue) = w0_sdplp25*vint_ci(intpos+1)+w1_sdplp25*vint_ci(intpos)
  intpos = intpos+3
end do
nlp_value = ilpvalue

end subroutine lp11_arblbr_ext_calcuvalue

subroutine lp12_arblbl_ext_calcuvalue(intentry,isma,nlp_value)

#include "drt_h.fh"
#include "pl_structure_h.fh"
#include "intsort_h.fh"
#include "lpextmode_h.fh"

w0_sdplp25 = (w0_sdplp-w1_sdplp)*w0g25
w1_sdplp25 = -2.d0*w0_sdplp*w0g25
intpos = intentry
ilpvalue = 0
do m_ia=1,nlsm_ext(isma)
  ilpvalue = ilpvalue+1
  value_lpext(ilpvalue) = w0_sdplp25*vint_ci(intpos+2)+w1_sdplp25*vint_ci(intpos)
  intpos = intpos+3

end do
nlp_value = ilpvalue

end subroutine lp12_arblbl_ext_calcuvalue

subroutine lp_arbr_ext_svtv_calcuvalue(intentry,nlp_value)

#include "drt_h.fh"
#include "lpextmode_h.fh"
#include "intsort_h.fh"

ivalue = 0
! G36a
w0lp = w0_plp*w0g36a
w1lp = w1_plp*w1g36a
valuelptmp1 = w0lp
w0lp = w0lp-w1lp
w1lp = valuelptmp1+w1lp
do intpos=intentry,intentry+ip2_intsymspace-3,3
  !intpos = intentry+lpext_int_index(ii)*3-3
  ivalue = ivalue+1
  ! ArBr -- B^rA^r =10
  value_lpext(ivalue) = vint_ci(intpos+2)*w0lp+vint_ci(intpos+1)*w1lp
end do
! G36b
! G1415
if (logic_g13) then
  w0lp = w0g13a*(w0_plp+w1_plp)
  intpos = intentry+ip2_intsymspace
  do ia=1,norb_ext
    ivalue = ivalue+1
    value_lpext(ivalue) = w0lp*vint_ci(intpos) !-vint_ci(intpos+1)*2.d0)
    intpos = intpos+2
  end do
end if
nlp_value = ivalue

end subroutine lp_arbr_ext_svtv_calcuvalue

subroutine lp_drr_ext_svtv_calcuvalue(intentry,nlp_value)

#include "drt_h.fh"
#include "lpextmode_h.fh"
#include "intsort_h.fh"

ivalue = 0
!G36a
w0lp = (w0_plp+w1_plp)*w0g36a
do intpos=intentry,intentry+ip2_drl_drl_intspace-2,2         !seve
  !intpos = intentry+lpext_int_index(ii)*3-3
  ivalue = ivalue+1
  value_lpext(ivalue) = vint_ci(intpos)*w0lp !-vint_ci(intpos)*w1lp      !severe_error_1111
end do

if (logic_g13) then
  !intpos = intentry+ip2_drl_drlintstart
  w0lp = 0.5d0*w0_plp*w0g13a
  !isma = ism_g1415
  !na_tmp = nlsm_ext(isma)
  do intpos=intentry+ip2_drl_drl_intspace,intentry+ip2_drl_intspace-2,2         !severe_new_e
    !do ia=1,na_tmp
    ivalue = ivalue+1
    value_lpext(ivalue) = vint_ci(intpos+1)*w0lp !+vint_ci(intpos+1)*w1lp
    !intpos = intpos+2
  end do
end if
nlp_value = ivalue

end subroutine lp_drr_ext_svtv_calcuvalue

subroutine lp9_drlbl_ext_calcuvalue_wyb(lri,lrk,isma)

#include "drt_h.fh"
#include "intsort_h.fh"
#include "lpextmode_h.fh"

next_sta = ibsm_ext(isma)-1
w0_sdplp25 = (w0_sdplp-w1_sdplp)*w0g25
w1_sdplp25 = -2.d0*w0_sdplp*w0g25
ia0 = (lri-1)*norb_ext
idorbint = lrk*2-2
ilpvalue = 0
do m_ia=1,nlsm_ext(isma)
  ilpvalue = ilpvalue+1
  ira = m_ia+next_sta
  ia = ia0+ira
  intposbase = intind_iaqq(ia)
  intpos = intposbase+idorbint
  value_lpext(ilpvalue) = w0_sdplp25*vint_ci(intpos)+w1_sdplp25*vint_ci(intpos+1)
end do

end subroutine lp9_drlbl_ext_calcuvalue_wyb

subroutine lp8_drlbr_sum_calcuvalue_wyb(lri,lrp,lrq,isma,nv)

#include "drt_h.fh"
#include "intsort_h.fh"
#include "pl_structure_h.fh"
#include "lpextmode_h.fh"
dimension vint_0(MAX_INNORB,MAX_EXTORB),vint_1(MAX_INNORB,MAX_EXTORB)

vint_0(1:norb_inn,1:norb_ext) = viasum_0(1:norb_inn,1:norb_ext)
vint_1(1:norb_inn,1:norb_ext) = viasum_1(1:norb_inn,1:norb_ext)

next_sta = ibsm_ext(isma)-1
ia0 = (lri-1)*norb_ext
idorbint_q = lrq*2-2
do m_ia=1,nlsm_ext(isma)
  ira = m_ia+next_sta
  ia = ia0+ira
  intposbase = intind_iaqq(ia)
  if (lrp /= 0) then
    idorbint_p = lrp*2-2
    intpos = intposbase+idorbint_p
    vint_0(lri,ira) = vint_0(lri,ira)-vint_ci(intpos)
    vint_1(lri,ira) = vint_1(lri,ira)-vint_ci(intpos+1)
  end if
  if (lrq /= 0) then
    idorbint_q = lrq*2-2
    intpos = intposbase+idorbint_q
    vint_0(lri,ira) = vint_0(lri,ira)-vint_ci(intpos)
    vint_1(lri,ira) = vint_1(lri,ira)-vint_ci(intpos+1)
  end if
end do

w0_sdplp25 = w0_sdplp*w0g25
w1_sdplp25 = 2*w0_sdplp*w0g25

ilpvalue = 0
do m_ia=1,nlsm_ext(isma)
  ilpvalue = ilpvalue+1
  ira = m_ia+next_sta
  value_lpext(ilpvalue) = w0_sdplp25*vint_0(lri,ira)-w1_sdplp25*vint_1(lri,ira) !value_lptm
end do
nv = ilpvalue

end subroutine lp8_drlbr_sum_calcuvalue_wyb

subroutine lp9_drlbl_sum_calcuvalue_wyb(lri,lrp,lrq,isma,nv)

#include "drt_h.fh"
#include "intsort_h.fh"
#include "lpextmode_h.fh"
dimension vint_0(MAX_INNORB,MAX_EXTORB),vint_1(MAX_INNORB,MAX_EXTORB)

vint_0(1:norb_INN,1:norb_ext) = viasum_0(1:norb_INN,1:norb_ext)
vint_1(1:norb_INN,1:norb_ext) = viasum_1(1:norb_INN,1:norb_ext)

next_sta = ibsm_ext(isma)-1
ia0 = (lri-1)*norb_ext
idorbint_q = lrq*2-2
do m_ia=1,nlsm_ext(isma)
  ira = m_ia+next_sta
  ia = ia0+ira
  intposbase = intind_iaqq(ia)
  if (lrp /= 0) then
    idorbint_p = lrp*2-2
    intpos = intposbase+idorbint_p
    vint_0(lri,ira) = vint_0(lri,ira)-vint_ci(intpos)
    vint_1(lri,ira) = vint_1(lri,ira)-vint_ci(intpos+1)
  end if
  if (lrq /= 0) then
    idorbint_q = lrq*2-2
    intpos = intposbase+idorbint_q
    vint_0(lri,ira) = vint_0(lri,ira)-vint_ci(intpos)
    vint_1(lri,ira) = vint_1(lri,ira)-vint_ci(intpos+1)
  end if
end do

next_sta = ibsm_ext(isma)-1
w0_sdplp25 = w0_sdplp*w0g25
w1_sdplp25 = -2.d0*w0_sdplp*w0g25

ilpvalue = 0
do m_ia=1,nlsm_ext(isma)
  ilpvalue = ilpvalue+1
  ira = m_ia+next_sta
  value_lpext(ilpvalue) = w0_sdplp25*vint_0(lri,ira)+w1_sdplp25*vint_1(lri,ira)
end do
nv = ilpvalue

end subroutine lp9_drlbl_sum_calcuvalue_wyb

subroutine lp_drl_ext_dd_calcuvalue_wyb(lri,iml,nlp_value)

#include "drt_h.fh"
#include "lpextmode_h.fh"
#include "intsort_h.fh"

w0lp = w0_plp*w0gdd
w1lp = w1_plp*w1gdd
nliml = nlsm_ext(iml)
iasta = ibsm_ext(iml)
iaend = iesm_ext(iml)
intpos0 = intind_abkk(lri)
ivalue = 0
jvalue = 0
if (logic_g49b) then
  do ira=iasta,iaend
    ivalue = ivalue+1
    lra = norb_number(ira)
    value_lpext(ivalue) = -w1lp*voint(lra,lri)
  end do
end if
mloop = nliml*(nliml-1)/2
intpos = intpos0+int_dd_drl*2
ivalue = ivalue+int_dd_drl
do I=1,mloop
  ivalue = ivalue+1
  jvalue = ivalue+mloop
  value_lpext(ivalue) = w0lp*vint_ci(intpos+1)-w1lp*vint_ci(intpos)
  value_lpext(jvalue) = value_lpext(ivalue)
  intpos = intpos+2
end do
nlp_value = jvalue

end subroutine lp_drl_ext_dd_calcuvalue_wyb

subroutine lp_arbl_ext_dd_calcuvalue(lri,lrj,iml,imr,nlp_value)

#include "drt_h.fh"
#include "lpextmode_h.fh"
#include "intsort_h.fh"

ij = LRI-norb_frz+NGW2(LRJ-norb_frz)
nliml = nlsm_ext(iml)
nlbf = ibsm_ext(iml)-1
intoffset = 2*nlbf
nlimr = nlsm_ext(imr)
w0lp = w0_plp*w0gdd
w1lp = w1_plp*w1gdd
valuetmp1 = w0lp
w0lp = w0lp-w1lp
w1lp = -valuetmp1*2.d0
ivalue = 0
! G50
if (logic_g50) then
  intentry = intind_ijcc(ij)

  mcloop = nliml
  intpos = intentry+intoffset
  do I=1,mcloop
    ivalue = ivalue+1
    value_lpext(ivalue) = vint_ci(intpos)*w0lp+vint_ci(intpos+1)*w1lp
    intpos = intpos+2
  end do
end if

intentry = intind_ijab(ij)
mcloop = nliml*nlimr
if (iml == imr) mcloop = nliml*(nliml-1)/2
INTPOS = INTENTRY
ivalue = ivalue
intoffset = int_dd_drl*3
INTPOS = INTENTRY+intoffset
ivalue = ivalue+int_dd_drl
if (logic_g49a) then
  ! G49a:Bl_Ar  line=12
  do I=1,MCLOOP
    ivalue = ivalue+1
    value_lpext(ivalue) = vint_ci(intpos+2)*w0lp+vint_ci(intpos)*w1lp
    INTPOS = INTPOS+3
  end do
end if
! G49b:Br_Al  line=11
INTPOS = INTENTRY+intoffset
if (logic_g49b) then
  do I=1,MCLOOP
    ivalue = ivalue+1
    value_lpext(ivalue) = vint_ci(intpos+1)*w0lp+vint_ci(intpos)*w1lp
    INTPOS = INTPOS+3
  end do
end if
nlp_value = ivalue

end subroutine lp_arbl_ext_dd_calcuvalue

subroutine lp_ar_coe_calcuvalue_wyb(idtu,isma,lri,lrj,nlp_value,lpcoe)

#include "drt_h.fh"
#include "intsort_h.fh"
#include "pl_structure_h.fh"
#include "lpextmode_h.fh"
real*8 valuetmp
dimension lpcoe(norb_dz+1:norb_inn)

! idtu=25,26,28,25(43),46,51,100(in act_space)

next_sta = ibsm_ext(isma)-1
w0_sdplp25 = w0_sdplp*w0g25
ia0 = (lri-1)*norb_ext
!idorbint = lri*2-2
!intpos = intposbase+idorbint
ilpvalue = 0
!intposbase = intentry

if (idtu == 100) then
  lsta = lri
  lend = norb_inn
  do m_ia=1,nlsm_ext(isma)
    ira = m_ia+next_sta
    lra = norb_number(ira)
    valuetmp = voint(lri,lra)
    !valuetmp = 0.D0
    ia = ia0+ira
    intposbase = intind_iaqq(ia)
    do iorb=lsta,lend
      KCOE = lpcoe(iorb)
      call NEOC(KCOE,NOCC,TCOE)
      idorbint = (iorb-1)*2
      intpos = intposbase+idorbint
      valuetmp = valuetmp+nocc*(vint_ci(intpos+1)+tcoe*vint_ci(intpos))
      !WL = WL+NEOC(K)*VLOP0*(VINT(LIST+2)+COE(K)*VINT(LIST+1))
    end do
    ilpvalue = ilpvalue+1
    value_lpext(ilpvalue) = w0_sdplp25*valuetmp
  end do
  nlp_value = nlsm_ext(isma)
  return
end if

if (idtu == 51) then
  ndorb = norb_DZ-lri
  do m_ia=1,nlsm_ext(isma)
    ira = m_ia+next_sta
    lra = norb_number(ira)
    valuetmp = voint(lri,lra)
    !valuetmp = 0.d0
    ia = ia0+ira
    intposbase = intind_iaqq(ia)
    idorbint = -2
    if (LOGIC_DH) then
      idorbint = lri*2-2
      intpos = intposbase+idorbint
      valuetmp = valuetmp+vint_ci(intpos)
      do idorb=1,ndorb
        idorbint = idorbint+2
        intpos = intposbase+idorbint
        valuetmp = valuetmp+2.d0*vint_ci(intpos+1)-vint_ci(intpos)
      end do
    end if
    iorbs = max(norb_dz+1,lri)
    do iorb=iorbs,norb_inn
      Kcoe = lpcoe(iorb)
      call NEOC(KCOE,NOCC,TCOE)
      idorbint = iorb*2-2
      intpos = intposbase+idorbint
      valuetmp = valuetmp+nocc*(vint_ci(intpos+1)+tcoe*vint_ci(intpos))
      !WL = WL+NEOC(K)*VLOP0*(VINT(LIST+2)+COE(K)*VINT(LIST+1))
    end do
    ilpvalue = ilpvalue+1
    value_lpext(ilpvalue) = w0_sdplp25*valuetmp
  end do
  nlp_value = nlsm_ext(isma)
  return
end if

if ((idtu == 25) .or. (idtu == 43)) then
  do m_ia=1,nlsm_ext(isma)
    ira = m_ia+next_sta
    lra = norb_number(ira)
    valuetmp = voint(lri,lra)
    !valuetmp = 0.d0
    ia = ia0+ira
    intposbase = intind_iaqq(ia)
    idorbint = lri*2-2
    intpos = intposbase+idorbint
    valuetmp = valuetmp+vint_ci(intpos)
    do idorb=lri+1,norb_dz
      idorbint = idorbint+2
      intpos = intposbase+idorbint
      valuetmp = valuetmp+2.d0*vint_ci(intpos+1)-vint_ci(intpos)
    end do
    iorbs = max(norb_dz+1,lri)
    do iorb=iorbs,norb_inn
      Kcoe = lpcoe(iorb)
      call NEOC(KCOE,NOCC,TCOE)
      idorbint = iorb*2-2
      intpos = intposbase+idorbint
      valuetmp = valuetmp+nocc*(vint_ci(intpos+1)+tcoe*vint_ci(intpos))
      !WL = WL+NEOC(K)*VLOP0*(VINT(LIST+2)+COE(K)*VINT(LIST+1))
    end do
    ilpvalue = ilpvalue+1
    value_lpext(ilpvalue) = w0_sdplp25*valuetmp
  end do
  nlp_value = nlsm_ext(isma)
  return
end if

if (idtu == 26) then
  ndorb = norb_DZ-lri
  do m_ia=1,nlsm_ext(isma)
    ira = m_ia+next_sta
    lra = norb_number(ira)
    valuetmp = voint(lri,lra)              !310
    ia = ia0+ira
    intposbase = intind_iaqq(ia)
    idorbint = lri*2-2
    if (LOGIC_DH) then
      do idorb=1,ndorb
        idorbint = idorbint+2
        intpos = intposbase+idorbint
        valuetmp = valuetmp+2.d0*vint_ci(intpos+1)-vint_ci(intpos)
      end do
    end if
    do iorb=norb_dz+1,norb_inn
      Kcoe = lpcoe(iorb)
      call NEOC(KCOE,NOCC,TCOE)
      idorbint = idorbint+2
      intpos = intposbase+idorbint
      valuetmp = valuetmp+nocc*(vint_ci(intpos+1)+tcoe*vint_ci(intpos))
      !WL = WL+NEOC(K)*VLOP0*(VINT(LIST+2)+COE(K)*VINT(LIST+1))
    end do

    ilpvalue = ilpvalue+1
    value_lpext(ilpvalue) = w0_sdplp25*valuetmp
  end do
  nlp_value = nlsm_ext(isma)
  return
end if
if ((idtu == 28) .or. (idtu == 46)) then
  nsorb = 1
  ndorb = norb_DZ-lri-1
  if (ndorb > 0) nsorb = 1
  if (idtu == 28) then
    icoe = -(JB_SYS+2)
  else !if (idtu == 46)
    icoe = 0
  end if
  do m_ia=1,nlsm_ext(isma)
    ira = m_ia+next_sta
    lra = norb_number(ira)
    ia = ia0+ira
    intposbase = intind_iaqq(ia)
    idorbint = lri*2-2
    intpos = intposbase+idorbint
    valuetmp = voint(lri,lra)+vint_ci(intpos)        ! 310+710
    if (nsorb > 0) then
      do iorb=lri+1,norb_dz
        idorbint = idorbint+2
        intpos = intposbase+idorbint
        if (iorb == lrj) valuetmp = valuetmp+vint_ci(intpos+1)+icoe*vint_ci(intpos)
        if (iorb /= lrj) valuetmp = valuetmp+2.d0*vint_ci(intpos+1)-vint_ci(intpos)
      end do
    end if
    do iorb=norb_dz+1,norb_inn
      Kcoe = lpcoe(iorb)
      call NEOC(KCOE,NOCC,TCOE)
      idorbint = idorbint+2
      intpos = intposbase+idorbint
      valuetmp = valuetmp+nocc*(vint_ci(intpos+1)+tcoe*vint_ci(intpos))
      !WL = WL+NEOC(K)*VLOP0*(VINT(LIST+2)+COE(K)*VINT(LIST+1))
    end do

    ilpvalue = ilpvalue+1
    value_lpext(ilpvalue) = w0_sdplp25*valuetmp
    intposbase = intposbase+norb_inn*2
  end do
  nlp_value = nlsm_ext(isma)
  return
end if
if ((idtu == 57) .or. (idtu == 29)) then
  nsorb = 1
  ndorb = norb_DZ-lri-1
  if (ndorb > 0) nsorb = 1
  if (idtu == 29) then
    icoe = JB_SYS
  else !if (idtu == 57)
    icoe = -1
  end if
  do m_ia=1,nlsm_ext(isma)
    ira = m_ia+next_sta
    lra = norb_number(ira)
    ia = ia0+ira
    intposbase = intind_iaqq(ia)
    idorbint = lri*2-2
    intpos = intposbase+idorbint
    valuetmp = voint(lri,lra)+vint_ci(intpos)        ! 310+710
    if (nsorb > 0) then
      do iorb=lri+1,norb_dz
        idorbint = idorbint+2
        intpos = intposbase+idorbint
        if (iorb == lrj) valuetmp = valuetmp+vint_ci(intpos+1)+icoe*vint_ci(intpos)
        if (iorb /= lrj) valuetmp = valuetmp+2.d0*vint_ci(intpos+1)-vint_ci(intpos)
      end do
    end if
    do iorb=norb_dz+1,norb_inn
      Kcoe = lpcoe(iorb)
      call NEOC(KCOE,NOCC,TCOE)
      idorbint = idorbint+2
      intpos = intposbase+idorbint
      valuetmp = valuetmp+nocc*(vint_ci(intpos+1)+tcoe*vint_ci(intpos))
      !WL = WL+NEOC(K)*VLOP0*(VINT(LIST+2)+COE(K)*VINT(LIST+1))
    end do

    ilpvalue = ilpvalue+1
    value_lpext(ilpvalue) = w0_sdplp25*valuetmp
    intposbase = intposbase+norb_inn*2
  end do
  nlp_value = nlsm_ext(isma)
  return
end if
!=======================================================================
if ((idtu == 55) .or. (idtu == 73)) then    !(11)(23)=55 (11)(13)=75
  do m_ia=1,nlsm_ext(isma)
    ira = m_ia+next_sta
    lra = norb_number(ira)
    valuetmp = voint(lri,lra)
    !valuetmp = 0.d0
    ia = ia0+ira
    intposbase = intind_iaqq(ia)
    idorbint = lri*2-2
    intpos = intposbase+idorbint
    valuetmp = valuetmp+vint_ci(intpos)
    do idorb=lri+1,norb_dz
      idorbint = idorbint+2
      intpos = intposbase+idorbint
      valuetmp = valuetmp+2.d0*vint_ci(intpos+1)-vint_ci(intpos)
    end do
    iorbs = max(norb_dz+1,lri)
    do iorb=iorbs,norb_inn
      Kcoe = lpcoe(iorb)
      call NEOC(KCOE,NOCC,TCOE)
      idorbint = iorb*2-2
      intpos = intposbase+idorbint
      valuetmp = valuetmp+nocc*(vint_ci(intpos+1)+tcoe*vint_ci(intpos))
      !WL = WL+NEOC(K)*VLOP0*(VINT(LIST+2)+COE(K)*VINT(LIST+1))
    end do
    ilpvalue = ilpvalue+1
    value_lpext(ilpvalue) = w0_sdplp25*valuetmp
  end do
  nlp_value = nlsm_ext(isma)
  return
end if

end subroutine lp_ar_coe_calcuvalue_wyb

subroutine lp_drl_ext_SS_calcuvalue(lri,nlp_value)

#include "drt_h.fh"
#include "lpextmode_h.fh"
#include "intsort_h.fh"

intpos = intind_abkk(lri)
intspace = intspace_abkk(lri)
ivalue = 0

! G2G4b
if (logic_g2g4a) then
  w0lp = w0_plp*w0g2a
  w1lp = w1_plp*w1g2a
  !ww0lp = w0_plp*w0g4a
  !ww1lp = w1_plp*w1g4a

  !valuelptmp1 = w0lp
  !w0lp = w0lp-w1lp
  !w1lp = -valuelptmp1*2.d0
  !valuelptmp1 = ww0lp
  !ww0lp = ww0lp-ww1lp
  !ww1lp = -valuelptmp1*2.d0

  do i=1,intspace
    ivalue = ivalue+2
    ! Drl -- B^lA^r =4_3
    value_lpext(ivalue) = vint_ci(intpos+1)*w0lp-vint_ci(intpos)*w1lp
    ! Drl -- B^rA^l =4_2
    value_lpext(ivalue-1) = (w0lp-w1lp)*vint_ci(intpos)
    intpos = intpos+2
  end do
end if

intpos = intind_abkk(lri)
w0lp = w0_plp*w0g36a
w1lp = w1_plp*w1g36a

do I=1,intspace
  ivalue = ivalue+1
  value_lpext(ivalue) = vint_ci(intpos+1)*w0lp-vint_ci(intpos)*w1lp
  intpos = intpos+2
end do
nlp_value = ivalue

end subroutine lp_drl_ext_SS_calcuvalue

subroutine lp_drl_sum_SS_calcuvalue(lri,lrj,nlp_value)

#include "drt_h.fh"
#include "lpextmode_h.fh"
#include "intsort_h.fh"
#include "pl_structure_h.fh"
dimension vint_0(norb_ext*norb_ext), vint_1(norb_ext*norb_ext)

intspace = intspace_abkk(1)
vint_0(1:intspace) = vijkk_0sum(1:intspace)
vint_1(1:intspace) = vijkk_1sum(1:intspace)

if (lri /= 0) then
  intpos = intind_abkk(lri)
  do I=1,intspace
    vint_0(I) = vint_0(I)-vint_ci(intpos)
    vint_1(I) = vint_1(I)-vint_ci(intpos+1)
    intpos = intpos+2
  end do
end if
if (lrj /= 0) then
  intpos = intind_abkk(lrj)
  do I=1,intspace
    vint_0(I) = vint_0(I)-vint_ci(intpos)
    vint_1(I) = vint_1(I)-vint_ci(intpos+1)
    intpos = intpos+2
  end do
end if

ivalue = 0
! G2G4b
if (logic_g2g4a) then
  w0lp = w0_plp*w0g2a
  w1lp = w1_plp*w1g2a
  do i=1,intspace
    ivalue = ivalue+2
    ! Drl -- B^lA^r =4_3
    value_lpext(ivalue) = vint_1(i)*w0lp-vint_0(i)*w1lp
    ! Drl -- B^rA^l =4_2
    value_lpext(ivalue-1) = (w0lp-w1lp)*vint_0(i)
  end do
end if

w0lp = w0_plp*w0g36a
w1lp = w1_plp*w1g36a

do I=1,intspace
  ivalue = ivalue+1
  value_lpext(ivalue) = vint_1(i)*w0lp-vint_0(i)*w1lp
end do

w0lp = w0_plp*w0g36a
w1lp = w1_plp*w1g36a

do I=1,intspace
  ivalue = ivalue+1
  value_lpext(ivalue) = vint_1(i)*w0lp-vint_0(i)*w1lp
end do
nlp_value = ivalue

end subroutine lp_drl_sum_SS_calcuvalue

subroutine lp_drl_ext_ST_calcuvalue(lri,nlp_value)

#include "drt_h.fh"
#include "lpextmode_h.fh"
#include "intsort_h.fh"

ivalue = 0
! G1415
if (logic_g1415) then

  w1lp = w1_plp*w1g14a

  do ismb=1,ng_sm
    isma = mul_tab(ismb,ism_g1415)
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
        !value_lpext(ivalue) = vint_ci(intposia)*ww0lp+vint_ci(intposia+1)*ww1lp+valuelpib
        ! OK,only for Spin=0
        value_lpext(ivalue) = (voint(lra,lri)-voint(lrb,lri))*w1lp
      end do
    end do
  end do
end if

intpos = intind_abkk(lri)
intspace = intspace_abkk(lri)
! G2G4b
if (logic_g2g4b) then
  w1lp = w1_plp*w1g4b

  do i=1,intspace
    ivalue = ivalue+2
    ! Drl -- B^lA^r =4_3
    value_lpext(ivalue) = w1lp*vint_ci(intpos)
    ! Drl -- B^rA^l =4_2
    value_lpext(ivalue-1) = -value_lpext(ivalue)
    intpos = intpos+2
  end do
end if

intpos = intind_abkk(lri)
w1lp = w1_plp*w1g36a

do I=1,intspace
  ivalue = ivalue+1
  value_lpext(ivalue) = -vint_ci(intpos)*w1lp
  intpos = intpos+2
end do
nlp_value = ivalue

end subroutine lp_drl_ext_ST_calcuvalue

subroutine lp_drl_ext_TT_calcuvalue(lri,n1415_value,nlp_value)

#include "drt_h.fh"
#include "lpextmode_h.fh"
#include "intsort_h.fh"

ivalue = 0
! G1415
if (logic_g1415) then

  w014 = w0_plp*w0g14a
  w114 = w1_plp*w1g14a
  w015 = w0_plp*w0g15a
  w115 = w1_plp*w1g15a

  w14lp = w014-w114
  w15lp = w015-w115

  do ismb=1,ng_sm
    isma = mul_tab(ismb,ism_g1415)
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
        value_lpext(ivalue) = w15lp*voint(lra,lri)+w14lp*voint(lrb,lri)
      end do
    end do
  end do
end if

n1415_value = ivalue
intpos = intind_abkk(lri)
intspace = intspace_abkk(lri)

w0lp = w0_plp*w0g36a
w1lp = w1_plp*w1g36a

do I=1,intspace
  ivalue = ivalue+1
  value_lpext(ivalue) = vint_ci(intpos+1)*w0lp-vint_ci(intpos)*w1lp
  intpos = intpos+2
end do
nlp_value = ivalue

end subroutine lp_drl_ext_TT_calcuvalue

subroutine lp_drl_sum_TT_calcuvalue(lri,lrj,n1415,nlp_value)

#include "drt_h.fh"
#include "lpextmode_h.fh"
#include "intsort_h.fh"
#include "pl_structure_h.fh"
dimension vint_0(norb_ext*norb_ext), vint_1(norb_ext*norb_ext)

intspace = intspace_abkk(1)
vint_0(1:intspace) = vijkk_0sum(1:intspace)
vint_1(1:intspace) = vijkk_1sum(1:intspace)

if (lri /= 0) then
  intpos = intind_abkk(lri)
  do I=1,intspace
    vint_0(I) = vint_0(I)-vint_ci(intpos)
    vint_1(I) = vint_1(I)-vint_ci(intpos+1)
    intpos = intpos+2
  end do
end if
if (lrj /= 0) then
  intpos = intind_abkk(lrj)
  do I=1,intspace
    vint_0(I) = vint_0(I)-vint_ci(intpos)
    vint_1(I) = vint_1(I)-vint_ci(intpos+1)
    intpos = intpos+2
  end do
end if

lrk0 = 1
if (lri == 1) lrk0 = 2
if (lrj == 2) lrk0 = 3
ivalue = 0
! G1415
if (logic_g1415) then

  w014 = w0_plp*w0g14a
  w015 = w0_plp*w0g15a

  w14lp = w014
  w15lp = w015

  do ismb=1,ng_sm
    isma = mul_tab(ismb,ism_g1415)
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
        value_lpext(ivalue) = w15lp*voint(lra,lrk0)+w14lp*voint(lrb,lrk0)
        do lrk=lrk0+1,norb_dz
          if (lrk == lri) cycle
          if (lrk == lrj) cycle
          value_lpext(ivalue) = value_lpext(ivalue)+w15lp*voint(lra,lrk)+w14lp*voint(lrb,lrk)
        end do
      end do
    end do
  end do
end if

n1415 = ivalue
w0lp = w0_plp*w0g36a

do I=1,intspace
  ivalue = ivalue+1
  value_lpext(ivalue) = vint_1(I)*w0lp
end do
nlp_value = ivalue

end subroutine lp_drl_sum_TT_calcuvalue

subroutine lp_arbr_ext_svtv_calcuvalue_wyb(LRI,LRJ,nlp_value)

#include "drt_h.fh"
#include "lpextmode_h.fh"
#include "intsort_h.fh"

IJ = LRI-NORB_FRZ+NGW2(LRJ-NORB_FRZ)
intentry = intind_ijab(ij)
intspace = intspace_ijab(ij)
ivalue = 0
intpos = intentry
! G36a
w0lp = w0_plp*w0g36a
w1lp = w1_plp*w1g36a
valuelptmp1 = w0lp
w0lp = w0lp-w1lp
w1lp = valuelptmp1+w1lp
! ArBr -- B^rA^r =10
do i=1,intspace
  ivalue = ivalue+1
  value_lpext(ivalue) = vint_ci(intpos+2)*w0lp+vint_ci(intpos+1)*w1lp
  intpos = intpos+3
end do
! G36b
! G1415
if (logic_g13) then
  intentry = intind_ijcc(IJ)
  intpos = intentry
  intspace = intspace_ijcc(ij)
  w0lp = w0g13a*(w0_plp+w1_plp)
  do i=1,intspace
    ivalue = ivalue+1
    ! ArBr -- D^r^r
    value_lpext(ivalue) = w0lp*vint_ci(intpos)
    intpos = intpos+2
  end do
end if
nlp_value = ivalue

end subroutine lp_arbr_ext_svtv_calcuvalue_wyb

subroutine lp_drr_ext_svtv_calcuvalue_wyb(lri,nlp_value)

#include "drt_h.fh"
#include "lpextmode_h.fh"
#include "intsort_h.fh"

intpos = intind_abkk(lri)
intspace = intspace_abkk(lri)
ivalue = 0
! G36a
w0lp = (w0_plp+w1_plp)*w0g36a
do I=1,intspace
  ivalue = ivalue+1
  value_lpext(ivalue) = vint_ci(intpos)*w0lp
  intpos = intpos+2
end do

if (logic_g13) then
  w0lp = 0.5d0*w0_plp*w0g13a
  do LRA=NORB_ALL,NORB_INN+1,-1
    ivalue = ivalue+1
    value_lpext(ivalue) = voint(LRA,LRI)*w0lp
  end do
end if
nlp_value = ivalue

end subroutine lp_drr_ext_svtv_calcuvalue_wyb

subroutine lp_arbl_ext_st_calcuvalue(lri,lrj,nlp_value)

#include "drt_h.fh"
#include "lpextmode_h.fh"
#include "intsort_h.fh"

ivalue = 0
ij = LRI-norb_frz+NGW2(LRJ-norb_frz)
intpos = intind_ijcc(IJ)
INTSPACE = intspace_ijcc(IJ)
!lmij = mul_tab(lsm_inn(lri),lsm_inn(lrj))
! G1415
if (logic_g1415) then

  w0lp = w0_plp*w0g14a
  w1lp = w1_plp*w1g14a
  ww0lp = w0_plp*w0g15a
  ww1lp = w1_plp*w1g15a

  valuelptmp1 = w0lp
  w0lp = w0lp-w1lp
  w1lp = -valuelptmp1*2.d0
  valuelptmp1 = ww0lp
  ww0lp = ww0lp-ww1lp
  ww1lp = -valuelptmp1*2.d0

  do ismb=1,ng_sm
    isma = mul_tab(ismb,ism_g1415)
    if (isma > ismb) cycle
    !lmab = mul_tab(isma,ismb)
    !if (lmab /= lmij) cycle
    ibsta = ibsm_ext(ismb)
    ibend = iesm_ext(ismb)
    iasta = ibsm_ext(isma)
    iaend = iesm_ext(isma)
    if (ismb == isma) ibsta = ibsta+1
    do ib=ibsta,ibend
      intposib = intpos+ib*2-2
      valuelpib = vint_ci(intposib)*w0lp+vint_ci(intposib+1)*w1lp
      do ia=iasta,min(iaend,ib-1)
        intposia = intpos+ia*2-2
        ivalue = ivalue+1
        valp = vint_ci(intposia)*ww0lp+vint_ci(intposia+1)*ww1lp
        value_lpext(ivalue) = valp+valuelpib
      end do
    end do
  end do
end if
! G13
if (logic_g13) then
  w0lp = w0g13a*w0_plp
  do ia=1,norb_ext
    intpos13 = intpos+ia*2-2
    ivalue = ivalue+1
    value_lpext(ivalue) = w0lp*(vint_ci(intpos13)-vint_ci(intpos13+1)*2.d0)
    intpos13 = intpos13+2
  end do
end if

INTSPACE = intspace_ijab(IJ)
! G2G4a
if (logic_g2g4a) then
  intpos = intind_ijab(IJ)

  w0lp = w0_plp*w0g2a
  w1lp = w1_plp*w1g2a
  ww0lp = w0_plp*w0g4a
  ww1lp = w1_plp*w1g4a

  valuelptmp1 = w0lp
  w0lp = w0lp-w1lp
  w1lp = -valuelptmp1*2.d0
  valuelptmp1 = ww0lp
  ww0lp = ww0lp-ww1lp
  ww1lp = -valuelptmp1*2.d0

  do i=1,intspace
    ivalue = ivalue+2
    ! ArBl -- B^lA^r =12
    value_lpext(ivalue-1) = vint_ci(intpos+2)*w0lp+vint_ci(intpos)*w1lp
    ! ArBl -- B^rA^l =11
    value_lpext(ivalue) = vint_ci(intpos+1)*ww0lp+vint_ci(intpos)*ww1lp
    intpos = intpos+3
  end do
else
  ! G2G4b
  if (logic_g2g4b) then
    intpos = intind_ijab(IJ)
    w0lp = w0_plp*w0g2b
    w1lp = w1_plp*w1g2b
    ww0lp = w0_plp*w0g4b
    ww1lp = w1_plp*w1g4b

    valuelptmp1 = w0lp
    w0lp = w0lp-w1lp
    w1lp = -valuelptmp1*2.d0
    valuelptmp1 = ww0lp
    ww0lp = ww0lp-ww1lp
    ww1lp = -valuelptmp1*2.d0

    do i=1,intspace
      ivalue = ivalue+2
      ! ArBl -- B^lA^r =12
      value_lpext(ivalue-1) = vint_ci(intpos+2)*ww0lp+vint_ci(intpos)*ww1lp
      ! ArBl -- B^rA^l =11
      value_lpext(ivalue) = vint_ci(intpos+1)*w0lp+vint_ci(intpos)*w1lp
      intpos = intpos+3
    end do
  end if

end if
! G36a
intpos = intind_ijab(IJ)
w0lp = w0_plp*w0g36a
w1lp = w1_plp*w1g36a
valuelptmp1 = w0lp
w0lp = w0lp-w1lp
w1lp = -valuelptmp1*2.d0
do i=1,intspace
  ivalue = ivalue+1
  ! ArBl -- B^lA^r =12
  value_lpext(ivalue) = vint_ci(intpos+2)*w0lp+vint_ci(intpos)*w1lp
  intpos = intpos+3
end do
! G36b
intpos = intind_ijab(IJ)
w0lp = w0_plp*w0g36b
w1lp = w1_plp*w1g36b
valuelptmp1 = w0lp
w0lp = w0lp-w1lp
w1lp = -valuelptmp1*2.d0
do i=1,intspace
  ivalue = ivalue+1
  ! ArBl -- B^rA^l =11
  value_lpext(ivalue) = vint_ci(intpos+1)*w0lp+vint_ci(intpos)*w1lp
  intpos = intpos+3
end do
nlp_value = ivalue

end subroutine lp_arbl_ext_st_calcuvalue

! lp7_ar_drl,lp7_drr_br,lp8_drl_br

subroutine lp678_ext_wyb_calcuvalue(lri,lrk,isma,nlp_value)

#include "drt_h.fh"
#include "intsort_h.fh"
#include "lpextmode_h.fh"

w0_sdplp25 = w0_sdplp*w0g25
ia0 = (lri-1)*norb_ext
intoffset = (lrk-1)*2
ilpvalue = 0
iasta = ibsm_ext(isma)
iaend = iesm_ext(isma)
do ia=iasta,iaend
  iaqq = ia0+ia
  intpos = intind_iaqq(iaqq)
  iposint = intpos+intoffset
  ilpvalue = ilpvalue+1
  value_lpext(ilpvalue) = w0_sdplp25*vint_ci(iposint)
end do
nlp_value = nlsm_ext(isma)

end subroutine lp678_ext_wyb_calcuvalue
