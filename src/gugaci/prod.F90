!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!                                                                      *
! Copyright (C) 2007, Bingbing Suo                                     *
!***********************************************************************
! 26 feb 2007 -bsuo- revised by suo bing for multi-root calculation

subroutine gdv_sequence_extspace(ilw,irw)
!***********************************************************************
! 26 feb 2007 - revised

use gugaci_global, only: ilsegdownwei, indx, log_prod, mcroot, value_lpext, vector1, vector2
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: ilw, irw
integer(kind=iwp) :: iij, irot, irtidx, mm, nn
real(kind=wp) :: valuelp, valuelptmp1, vlptmp

!write(u6,*) '  dv_test ','  vd_test '

! mrpt2
if (log_prod == 3) then
  call gdv_sequence_extspace_pt(ilw,irw)
  return
end if

do irot=1,mcroot
  irtidx = indx(irot)
  mm = ilw+irtidx
  nn = irw+1+irtidx
  valuelp = vector2(nn)
  valuelptmp1 = vector1(nn)
  do iij=1,ilsegdownwei
    vlptmp = value_lpext(iij)
    mm = mm+1
    vector2(mm) = vector2(mm)+valuelptmp1*vlptmp
    valuelp = valuelp+vector1(mm)*vlptmp
    !if ((mm == 86) .and. (nn == 5)) then
    !  write(u6,*) 'dv_test 1 ',vlptmp,vector2(mm)
    !end if
    !if ((mm == 5) .and. (nn == 86)) then
    !  write(u6,*) 'dv_test 2 ',valuelp,vector2(nn)
    !end if
  end do
  vector2(nn) = valuelp
end do

end subroutine gdv_sequence_extspace

subroutine gdv_sequence_extspace_pt(ilw,irw)  !log_prod=

use gugaci_global, only: ilsegdownwei, value_lpext, vcm, vector2
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: ilw, irw
integer(kind=iwp) :: iij, mm, nn
real(kind=wp) :: wl

mm = ilw
nn = irw+1
do iij=1,ilsegdownwei
  mm = mm+1
  wl = value_lpext(iij)
  vector2(mm) = vector2(mm)+vcm(nn)*wl
  !if ((jpad == 1) .and. (ipae == 2)) then
  !  write(u6,'(a,4i6,3f18.8)') 'dv_test_pt 1 ',jpadl,ipael,mm,nn,wl,vcm(nn),vector2(mm)
  !end if
  !if ((jpadl == 1) .and. (ipael == 2)) then
  !  write(u6,'(a,4i6,3f18.8)') 'dv_test_pt 2 ',jpad,ipae,mm,nn,wl,vcm(nn),vector2(mm)
  !end if
end do

end subroutine gdv_sequence_extspace_pt

subroutine gtd_sequence_extspace(iplplweiin,iplprweiin)
!***********************************************************************
! 26 feb 2007 - revised

use gugaci_global, only: indx, iweista_g25, iweista_g28, logic_g25a, logic_g25b, logic_g28a, mcroot, nint_g25, nint_g28, nwei_g25, &
                         nwei_g28, value_lpext, vector1, vector2
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: iplplweiin, iplprweiin
integer(kind=iwp) :: i, ilpvalue, iplplwei, iplprwei, irot, irtidx, itmp, mm, nn, nn0
real(kind=wp) :: vlptmp, vlptmp1, vlptmp2

!write(u6,*) ' td_test _1/2',' dt_test '

do irot=1,mcroot
  irtidx = indx(irot)
  iplplwei = iplplweiin+irtidx
  iplprwei = iplprweiin+irtidx

  ilpvalue = 0
  if (logic_g25a) then
    mm = iplplwei+iweista_g25-1
    nn0 = iplprwei
    do itmp=1,nint_g25
      ilpvalue = ilpvalue+1
      vlptmp = value_lpext(ilpvalue)
      nn = nn0
      do i=1,nwei_g25
        mm = mm+1
        nn = nn+1
        vector2(mm) = vector2(mm)+vector1(nn)*vlptmp
        vector2(nn) = vector2(nn)+vector1(mm)*vlptmp
      end do
    end do
  else if (logic_g25b) then
    mm = iplplwei+iweista_g25-1
    nn0 = iplprwei
    ilpvalue = ilpvalue+1
    do itmp=2,nint_g25
      ilpvalue = ilpvalue+1
      vlptmp = value_lpext(ilpvalue)
      nn = nn0
      do i=1,itmp-1
        mm = mm+1
        nn = nn+1
        vector2(mm) = vector2(mm)+vector1(nn)*vlptmp
        vector2(nn) = vector2(nn)+vector1(mm)*vlptmp
      end do
    end do
    mm = iplplwei+iweista_g28-1
    nn = iplprwei
    nn = nn+1
    do itmp=2,nwei_g28
      nn = nn+1
      vlptmp = vector2(nn)
      vlptmp1 = vector1(nn)
      ilpvalue = 0
      do i=1,itmp-1
        ilpvalue = ilpvalue+1
        vlptmp2 = -value_lpext(ilpvalue)
        mm = mm+1
        vector2(mm) = vector2(mm)+vlptmp1*vlptmp2
        vlptmp = vlptmp+vector1(mm)*vlptmp2
      end do
      vector2(nn) = vlptmp
    end do
  else if (logic_g28a) then
    mm = iplplwei+iweista_g28-1
    nn0 = iplprwei
    do nn=nn0+1,nn0+nwei_g28
      vlptmp = vector2(nn)
      vlptmp1 = vector1(nn)
      ilpvalue = 0
      do i=1,nint_g28
        ilpvalue = ilpvalue+1
        vlptmp2 = -value_lpext(ilpvalue)
        mm = mm+1
        vector2(mm) = vector2(mm)+vlptmp1*vlptmp2
        vlptmp = vlptmp+vector1(mm)*vlptmp2
      end do
      vector2(nn) = vlptmp
    end do
  end if
end do

end subroutine gtd_sequence_extspace

subroutine inn_ext_ss_loop_unpack(iplplweiin,iplprweiin)
!***********************************************************************
! 26 feb 2007 - revised

use gugaci_global, only: ibsm_ext, idownwei_g131415, iesm_ext, indx, ism_g2g4, iwt_sm_s_ext, logic_g1415, logic_g2g4a, &
                         logic_g2g4b, logic_g34a, logic_g34b, logic_g35a, logic_g35b, logic_g36a, logic_g36b, lpend34a, lpend34b, &
                         lpend35a, lpend35b, lpend36a, lpend36b, lpext_wei, lpsta34a, lpsta34b, lpsta35a, lpsta35b, lpsta36a, &
                         lpsta36b, mcroot, ng_sm, nvalue_space_ss, value_lpext, vector1, vector2
use Symmetry_Info, only: Mul
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: iplplweiin, iplprweiin
integer(kind=iwp) :: i, ia, iaend, iasta, ib, ibend, ibsta, icle, ii, ii0, iii, iij, ilwtmp, iplplwei, iplprwei, irot, irtidx, &
                     irwtmp, isma, ismb, lpend34, lpend35, lpend36, lpsta34, lpsta35, lpsta36, mm, mm0, nn, nn0, nna, nnb
real(kind=wp) :: valuelp, valuelp1, valuelp2, valuelptmp1, valuelptmp2, valuetmp

!write(u6,*) ' ss_test 1/2'
do irot=1,mcroot
  irtidx = indx(irot)
  iplplwei = iplplweiin+irtidx
  iplprwei = iplprweiin+irtidx

  ii = 1
  if (logic_g1415) then
    mm = iplplwei
    nn = iplprwei
    do i=1,idownwei_g131415
      valuelp = value_lpext(ii)
      mm = mm+1
      nn = nn+1
      vector2(mm) = vector2(mm)+vector1(nn)*valuelp
      vector2(nn) = vector2(nn)+vector1(mm)*valuelp
      ii = ii+1
    end do
  end if

  ii0 = ii
  if (logic_g2g4a) then
    ii = ii0
    mm0 = iplplwei
    nn0 = iplprwei+iwt_sm_s_ext
    mm = mm0
    do ismb=1,ng_sm
      isma = Mul(ismb,ism_g2g4)
      if (isma > ismb) cycle
      ibsta = ibsm_ext(ismb)
      ibend = iesm_ext(ismb)
      iasta = ibsm_ext(isma)
      iaend = iesm_ext(isma)
      if (ismb == isma) ibsta = ibsta+1
      nna = nn0+ibsta-1
      do ib=ibsta,ibend
        nna = nna+1
        nnb = nn0+iasta-1
        valuelptmp1 = vector1(nna)
        valuelptmp2 = vector2(nna)
        do ia=iasta,min(iaend,ib-1)
          valuelp1 = value_lpext(ii)
          valuelp2 = value_lpext(ii+1)
          mm = mm+1
          nnb = nnb+1
          vector2(mm) = vector2(mm)+valuelptmp1*valuelp1+vector1(nnb)*valuelp2
          vector2(nnb) = vector2(nnb)+vector1(mm)*valuelp2
          valuelptmp2 = valuelptmp2+vector1(mm)*valuelp1
          ii = ii+2
        end do
        vector2(nna) = valuelptmp2
      end do
    end do
  end if
  if (logic_g2g4b) then
    ii = ii0
    mm0 = iplprwei
    nn0 = iplplwei+iwt_sm_s_ext
    mm = mm0
    do ismb=1,ng_sm
      isma = Mul(ismb,ism_g2g4)
      if (isma > ismb) cycle
      ibsta = ibsm_ext(ismb)
      ibend = iesm_ext(ismb)
      iasta = ibsm_ext(isma)
      iaend = iesm_ext(isma)
      if (ismb == isma) ibsta = ibsta+1
      nna = nn0+ibsta-1
      do ib=ibsta,ibend
        nna = nna+1
        nnb = nn0+iasta-1
        valuelptmp1 = vector1(nna)
        valuelptmp2 = vector2(nna)
        do ia=iasta,min(iaend,ib-1)
          valuelp2 = value_lpext(ii)
          valuelp1 = value_lpext(ii+1)
          mm = mm+1
          nnb = nnb+1
          vector2(mm) = vector2(mm)+valuelptmp1*valuelp1+vector1(nnb)*valuelp2
          vector2(nnb) = vector2(nnb)+vector1(mm)*valuelp2
          valuelptmp2 = valuelptmp2+vector1(mm)*valuelp1
          ii = ii+2
        end do
        vector2(nna) = valuelptmp2
      end do
    end do
  end if
  !iaddii = (ii-ii0)/2
  ii0 = ii-1

  do icle=1,2
    if (((icle == 1) .and. logic_g36a) .or. ((icle == 2) .and. logic_g36b)) then
      if (icle == 1) then
        mm0 = iplplwei
        nn0 = iplprwei
        lpsta36 = lpsta36a
        lpend36 = lpend36a
      else
        mm0 = iplprwei
        nn0 = iplplwei
        lpsta36 = lpsta36b
        lpend36 = lpend36b
      end if
      do iii=lpsta36,lpend36,4
        ilwtmp = lpext_wei(iii)
        irwtmp = lpext_wei(iii+1)
        ii = ii0+lpext_wei(iii+2)
        mm = mm0+ilwtmp
        nn = nn0+irwtmp
        do i=1,lpext_wei(iii+3)
          valuetmp = value_lpext(ii)
          vector2(mm) = vector2(mm)+vector1(nn)*valuetmp
          vector2(nn) = vector2(nn)+vector1(mm)*valuetmp
          mm = mm+1
          nn = nn+1
        end do
      end do
    end if

    if (((icle == 1) .and. logic_g35a) .or. ((icle == 2) .and. logic_g35b)) then
      if (icle == 1) then
        mm0 = iplplwei
        nn0 = iplprwei
        lpsta35 = lpsta35a
        lpend35 = lpend35a
      else
        mm0 = iplprwei
        nn0 = iplplwei
        lpsta35 = lpsta35b
        lpend35 = lpend35b
      end if
      do iii=lpsta35,lpend35,4
        mm = mm0+lpext_wei(iii)
        nn = nn0+lpext_wei(iii+1)
        iij = ii0+lpext_wei(iii+2)
        valuelp = vector2(nn)
        valuelptmp1 = vector1(nn)
        do i=1,lpext_wei(iii+3)
          valuetmp = value_lpext(iij)
          vector2(mm) = vector2(mm)+valuelptmp1*valuetmp
          valuelp = valuelp+vector1(mm)*valuetmp
          iij = iij+1
          mm = mm+1
        end do
        vector2(nn) = valuelp
      end do
    end if

    !cycle

    if (((icle == 1) .and. logic_g34a) .or. ((icle == 2) .and. logic_g34b)) then
      if (icle == 1) then
        mm0 = iplplwei
        nn0 = iplprwei
        lpsta34 = lpsta34a
        lpend34 = lpend34a
      else
        mm0 = iplprwei
        nn0 = iplplwei
        lpsta34 = lpsta34b
        lpend34 = lpend34b
      end if
      do iii=lpsta34,lpend34,4
        ilwtmp = lpext_wei(iii)
        irwtmp = lpext_wei(iii+1)
        mm = mm0+ilwtmp
        nn = nn0+irwtmp
        iij = ii0+lpext_wei(iii+2)
        valuelp = vector2(nn)
        valuelptmp1 = vector1(nn)
        do i=1,lpext_wei(iii+3)
          valuetmp = value_lpext(iij)
          vector2(mm) = vector2(mm)+valuelptmp1*valuetmp
          valuelp = valuelp+vector1(mm)*valuetmp
          iij = iij+1
          mm = mm+1
        end do
        vector2(nn) = valuelp
      end do
    end if

    ii0 = ii0+nvalue_space_ss
  end do

end do

end subroutine inn_ext_ss_loop_unpack

subroutine inn_ext_ss_drl_loop_unpack(iplplweiin,iplprweiin)
!***********************************************************************
! 26 feb 2007 - revised

use gugaci_global, only: ibsm_ext, idownwei_g131415, iesm_ext, indx, ism_g2g4, iwt_sm_s_ext, logic_g1415, logic_g2g4a, &
                         logic_g2g4b, logic_g34a, logic_g34b, logic_g35a, logic_g35b, logic_g36a, logic_g36b, lpend34a, lpend34b, &
                         lpend35a, lpend35b, lpend36a, lpend36b, lpext_wei, lpsta34a, lpsta34b, lpsta35a, lpsta35b, lpsta36a, &
                         lpsta36b, mcroot, ng_sm, value_lpext, vector1, vector2
use Symmetry_Info, only: Mul
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: iplplweiin, iplprweiin
integer(kind=iwp) :: i, ia, iaend, iasta, ib, ibend, ibsta, icle, ii, ii0, iii, iij, ilwtmp, iplplwei, iplprwei, irot, irtidx, &
                     irwtmp, isma, ismb, lpend34, lpend35, lpend36, lpsta34, lpsta35, lpsta36, mm, mm0, nn, nn0, nna, nnb
real(kind=wp) :: valuelp, valuelp1, valuelp2, valuelptmp1, valuelptmp2, valuetmp

!write(u6,*) ' ss_test 2/2'
do irot=1,mcroot
  irtidx = indx(irot)
  iplplwei = iplplweiin+irtidx
  iplprwei = iplprweiin+irtidx

  ii = 1
  if (logic_g1415) then
    mm = iplplwei
    nn = iplprwei
    do i=1,idownwei_g131415
      valuelp = value_lpext(ii)
      mm = mm+1
      nn = nn+1
      vector2(mm) = vector2(mm)+vector1(nn)*valuelp
      vector2(nn) = vector2(nn)+vector1(mm)*valuelp
      ii = ii+1
    end do
  end if

  ii0 = ii
  if (logic_g2g4a) then
    ii = ii0
    mm0 = iplplwei
    nn0 = iplprwei+iwt_sm_s_ext
    mm = mm0
    do ismb=1,ng_sm
      isma = Mul(ismb,ism_g2g4)
      if (isma > ismb) cycle
      ibsta = ibsm_ext(ismb)
      ibend = iesm_ext(ismb)
      iasta = ibsm_ext(isma)
      iaend = iesm_ext(isma)
      if (ismb == isma) ibsta = ibsta+1
      nna = nn0+ibsta-1
      do ib=ibsta,ibend
        nna = nna+1
        nnb = nn0+iasta-1
        valuelptmp1 = vector1(nna)
        valuelptmp2 = vector2(nna)
        do ia=iasta,min(iaend,ib-1)
          valuelp1 = value_lpext(ii)
          valuelp2 = value_lpext(ii+1)
          mm = mm+1
          nnb = nnb+1
          vector2(mm) = vector2(mm)+valuelptmp1*valuelp2+vector1(nnb)*valuelp2
          vector2(nnb) = vector2(nnb)+vector1(mm)*valuelp2
          valuelptmp2 = valuelptmp2+vector1(mm)*valuelp2
          ii = ii+2
        end do
        vector2(nna) = valuelptmp2
      end do
    end do
  end if

  if (logic_g2g4b) then
    ii = ii0
    mm0 = iplprwei
    nn0 = iplplwei+iwt_sm_s_ext
    mm = mm0
    do ismb=1,ng_sm
      isma = Mul(ismb,ism_g2g4)
      if (isma > ismb) cycle
      ibsta = ibsm_ext(ismb)
      ibend = iesm_ext(ismb)
      iasta = ibsm_ext(isma)
      iaend = iesm_ext(isma)
      if (ismb == isma) ibsta = ibsta+1
      nna = nn0+ibsta-1
      do ib=ibsta,ibend
        nna = nna+1
        nnb = nn0+iasta-1
        valuelptmp1 = vector1(nna)
        valuelptmp2 = vector2(nna)
        do ia=iasta,min(iaend,ib-1)
          valuelp2 = value_lpext(ii)
          valuelp1 = value_lpext(ii+1)
          mm = mm+1
          nnb = nnb+1
          vector2(mm) = vector2(mm)+valuelptmp1*valuelp1+vector1(nnb)*valuelp1
          vector2(nnb) = vector2(nnb)+vector1(mm)*valuelp1
          valuelptmp2 = valuelptmp2+vector1(mm)*valuelp1
          ii = ii+2
        end do
        vector2(nna) = valuelptmp2
      end do
    end do
  end if
  !iaddii = (ii-ii0)/2
  ii0 = ii-1

  do icle=1,2
    if (((icle == 1) .and. logic_g36a) .or. ((icle == 2) .and. logic_g36b)) then
      if (icle == 1) then
        mm0 = iplplwei
        nn0 = iplprwei
        lpsta36 = lpsta36a
        lpend36 = lpend36a
      else
        mm0 = iplprwei
        nn0 = iplplwei
        lpsta36 = lpsta36b
        lpend36 = lpend36b
      end if
      do iii=lpsta36,lpend36,4
        ilwtmp = lpext_wei(iii)
        irwtmp = lpext_wei(iii+1)
        ii = ii0+lpext_wei(iii+2)
        mm = mm0+ilwtmp
        nn = nn0+irwtmp
        do i=1,lpext_wei(iii+3)
          valuetmp = value_lpext(ii)
          vector2(mm) = vector2(mm)+vector1(nn)*valuetmp
          vector2(nn) = vector2(nn)+vector1(mm)*valuetmp
          mm = mm+1
          nn = nn+1
        end do
      end do
    end if

    if (((icle == 1) .and. logic_g35a) .or. ((icle == 2) .and. logic_g35b)) then
      if (icle == 1) then
        mm0 = iplplwei
        nn0 = iplprwei
        lpsta35 = lpsta35a
        lpend35 = lpend35a
      else
        mm0 = iplprwei
        nn0 = iplplwei
        lpsta35 = lpsta35b
        lpend35 = lpend35b
      end if
      do iii=lpsta35,lpend35,4
        mm = mm0+lpext_wei(iii)
        nn = nn0+lpext_wei(iii+1)
        iij = ii0+lpext_wei(iii+2)
        valuelp = vector2(nn)
        valuelptmp1 = vector1(nn)
        do i=1,lpext_wei(iii+3)
          valuetmp = value_lpext(iij)
          vector2(mm) = vector2(mm)+valuelptmp1*valuetmp
          valuelp = valuelp+vector1(mm)*valuetmp
          iij = iij+1
          mm = mm+1
        end do
        vector2(nn) = valuelp
      end do
    end if

    !cycle

    if (((icle == 1) .and. logic_g34a) .or. ((icle == 2) .and. logic_g34b)) then
      if (icle == 1) then
        mm0 = iplplwei
        nn0 = iplprwei
        lpsta34 = lpsta34a
        lpend34 = lpend34a
      else
        mm0 = iplprwei
        nn0 = iplplwei
        lpsta34 = lpsta34b
        lpend34 = lpend34b
      end if
      do iii=lpsta34,lpend34,4
        ilwtmp = lpext_wei(iii)
        irwtmp = lpext_wei(iii+1)
        mm = mm0+ilwtmp
        nn = nn0+irwtmp
        iij = ii0+lpext_wei(iii+2)
        valuelp = vector2(nn)
        valuelptmp1 = vector1(nn)
        do i=1,lpext_wei(iii+3)
          valuetmp = value_lpext(iij)
          vector2(mm) = vector2(mm)+valuelptmp1*valuetmp
          valuelp = valuelp+vector1(mm)*valuetmp
          iij = iij+1
          mm = mm+1
        end do
        vector2(nn) = valuelp
      end do
    end if

  end do

end do

end subroutine inn_ext_ss_drl_loop_unpack

subroutine inn_ext_sv_loop_unpack(ilw,irw)   !,ilsegdownwei)
!***********************************************************************
! 26 feb 2007 - revised

use gugaci_global, only: ilsegdownwei, indx, log_prod, mcroot, value_lpext, vector1, vector2
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: ilw, irw
integer(kind=iwp) :: iij, irot, irtidx, mm, nn
real(kind=wp) :: valuelp, valuelptmp1

!write(u6,*) 'sv_test,   tv_test '
if (log_prod == 3) then
  call inn_ext_svloop_unpack_pt(ilw,irw)
  return
end if

do irot=1,mcroot
  irtidx = indx(irot)
  mm = ilw+irtidx
  nn = irw+1+irtidx
  valuelp = vector2(nn)
  valuelptmp1 = vector1(nn)
  do iij=1,ilsegdownwei
    mm = mm+1
    vector2(mm) = vector2(mm)+valuelptmp1*value_lpext(iij)
    valuelp = valuelp+vector1(mm)*value_lpext(iij)
  end do
  vector2(nn) = valuelp
end do

end subroutine inn_ext_sv_loop_unpack

subroutine inn_ext_svloop_unpack_pt(ilw,irw)

use gugaci_global, only: ilsegdownwei, value_lpext, vcm, vector2
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: ilw, irw
integer(kind=iwp) :: iij, mm, nn
real(kind=wp) :: wl
!character(len=16) :: loop_type

!loop_type = ' sv_test_pt, tv_test_pt'
mm = ilw
nn = irw+1
do iij=1,ilsegdownwei
  mm = mm+1
  wl = value_lpext(iij)
  vector2(mm) = vector2(mm)+vcm(nn)*wl
  !vector2(mm) = vector2(mm)+vcm0(nn)*wl
  !if (mm == mtest)then
  !  write(nf2,'(2i5,3f18.10)') 1, nn,vcm0(nn),wl,vector2(mm-ndim_h0)
  !end if
end do

end subroutine inn_ext_svloop_unpack_pt

subroutine inn_ext_tt_loop_unpack(iplplweiin,iplprweiin)
!***********************************************************************
! 26 feb 2007 - revised

use gugaci_global, only: idownwei_g131415, indx, logic_g1415, logic_g34a, logic_g34b, logic_g35a, logic_g35b, logic_g36a, &
                         logic_g36b, lpend34a, lpend34b, lpend35a, lpend35b, lpend36a, lpend36b, lpext_wei, lpsta34a, lpsta34b, &
                         lpsta35a, lpsta35b, lpsta36a, lpsta36b, mcroot, nvalue_space_ss, value_lpext, vector1, vector2
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: iplplweiin, iplprweiin
integer(kind=iwp) :: i, icle, ii, ii0, iii, iij, ilwtmp, iplplwei, iplprwei, irot, irtidx, irwtmp, lpend34, lpend35, lpend36, &
                     lpsta34, lpsta35, lpsta36, mm, mm0, nn, nn0
real(kind=wp) :: valuelp, valuelptmp1, valuetmp

!write(u6,*) '  tt_test 1/2'

do irot=1,mcroot
  irtidx = indx(irot)
  iplplwei = iplplweiin+irtidx
  iplprwei = iplprweiin+irtidx

  ii = 1
  if (logic_g1415) then
    mm = iplplwei
    nn = iplprwei
    do i=1,idownwei_g131415
      valuelp = value_lpext(ii)
      mm = mm+1
      nn = nn+1
      vector2(mm) = vector2(mm)+vector1(nn)*valuelp
      vector2(nn) = vector2(nn)+vector1(mm)*valuelp
      ii = ii+1
    end do
  end if

  ii0 = ii-1         !severe_error

  do icle=1,2
    if (((icle == 1) .and. logic_g36a) .or. ((icle == 2) .and. logic_g36b)) then
      if (icle == 1) then
        mm0 = iplplwei
        nn0 = iplprwei
        lpsta36 = lpsta36a
        lpend36 = lpend36a
      else
        mm0 = iplprwei
        nn0 = iplplwei
        lpsta36 = lpsta36b
        lpend36 = lpend36b
      end if
      do iii=lpsta36,lpend36,4
        ilwtmp = lpext_wei(iii)
        irwtmp = lpext_wei(iii+1)
        ii = ii0+lpext_wei(iii+2)
        mm = mm0+ilwtmp
        nn = nn0+irwtmp
        do i=1,lpext_wei(iii+3)
          valuelp = value_lpext(ii)
          vector2(mm) = vector2(mm)+vector1(nn)*valuelp
          vector2(nn) = vector2(nn)+vector1(mm)*valuelp
          mm = mm+1
          nn = nn+1
        end do
      end do
    end if

    if (((icle == 1) .and. logic_g35a) .or. ((icle == 2) .and. logic_g35b)) then
      if (icle == 1) then
        mm0 = iplplwei
        nn0 = iplprwei
        lpsta35 = lpsta35a
        lpend35 = lpend35a
      else
        mm0 = iplprwei
        nn0 = iplplwei
        lpsta35 = lpsta35b
        lpend35 = lpend35b
      end if
      do iii=lpsta35,lpend35,4
        mm = mm0+lpext_wei(iii)
        nn = nn0+lpext_wei(iii+1)
        iij = ii0+lpext_wei(iii+2)
        valuelp = vector2(nn)
        valuelptmp1 = vector1(nn)
        do i=1,lpext_wei(iii+3)
          valuetmp = -value_lpext(iij)
          vector2(mm) = vector2(mm)+valuelptmp1*valuetmp
          valuelp = valuelp+vector1(mm)*valuetmp
          iij = iij+1
          mm = mm+1
        end do
        vector2(nn) = valuelp      !severe_error_1202
      end do
    end if

    if (((icle == 1) .and. logic_g34a) .or. ((icle == 2) .and. logic_g34b)) then
      if (icle == 1) then
        mm0 = iplplwei
        nn0 = iplprwei
        lpsta34 = lpsta34a
        lpend34 = lpend34a
      else
        mm0 = iplprwei
        nn0 = iplplwei
        lpsta34 = lpsta34b
        lpend34 = lpend34b
      end if
      do iii=lpsta34,lpend34,4
        ilwtmp = lpext_wei(iii)
        irwtmp = lpext_wei(iii+1)
        mm = mm0+ilwtmp
        nn = nn0+irwtmp
        iij = ii0+lpext_wei(iii+2)
        valuelp = vector2(nn)
        valuelptmp1 = vector1(nn)
        do i=1,lpext_wei(iii+3)
          vector2(mm) = vector2(mm)+valuelptmp1*value_lpext(iij)
          valuelp = valuelp+vector1(mm)*value_lpext(iij)
          iij = iij+1
          mm = mm+1
        end do
        vector2(nn) = valuelp
      end do
    end if
    ii0 = ii0+nvalue_space_ss
  end do

end do

end subroutine inn_ext_tt_loop_unpack

subroutine inn_ext_ts_loop_unpack(iplplweiin,iplprweiin)
!***********************************************************************
! 26 feb 2007 - revised

use gugaci_global, only: ibsm_ext, idownwei_g131415, iesm_ext, indx, ism_g2g4, iwt_sm_s_ext, logic_g1415, logic_g2g4a, &
                         logic_g2g4b, logic_g34a, logic_g34b, logic_g35a, logic_g35b, logic_g36a, logic_g36b, lpend34a, lpend34b, &
                         lpend35a, lpend35b, lpend36a, lpend36b, lpext_wei, lpsta34a, lpsta34b, lpsta35a, lpsta35b, lpsta36a, &
                         lpsta36b, mcroot, ng_sm, nvalue_space_ss, value_lpext, vector1, vector2
use Symmetry_Info, only: Mul
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: iplplweiin, iplprweiin
integer(kind=iwp) :: i, ia, iaend, iasta, ib, ibend, ibsta, icle, ii, ii0, iii, iij, ilwtmp, iplplwei, iplprwei, irot, irtidx, &
                     irwtmp, isma, ismb, lpend34, lpend35, lpend36, lpsta34, lpsta35, lpsta36, mm, mm0, nn, nn0, nna, nnb
real(kind=wp) :: valuelp, valuelp1, valuelp2, valuelptmp1, valuelptmp2, valuetmp

!write(u6,*) '  ts_test 1/2 '

do irot=1,mcroot
  irtidx = indx(irot)
  iplplwei = iplplweiin+irtidx
  iplprwei = iplprweiin+irtidx

  ii = 1
  if (logic_g1415) then
    mm = iplplwei
    nn = iplprwei
    do i=1,idownwei_g131415
      valuelp = value_lpext(ii)
      mm = mm+1
      nn = nn+1
      vector2(mm) = vector2(mm)+vector1(nn)*valuelp
      vector2(nn) = vector2(nn)+vector1(mm)*valuelp
      ii = ii+1
    end do
  end if

  ii0 = ii
  if (logic_g2g4a) then
    ii = ii0
    mm0 = iplplwei
    nn0 = iplprwei+iwt_sm_s_ext
    mm = mm0            !severe_error
    do ismb=1,ng_sm
      isma = Mul(ismb,ism_g2g4)
      if (isma > ismb) cycle
      ibsta = ibsm_ext(ismb)
      ibend = iesm_ext(ismb)
      iasta = ibsm_ext(isma)
      iaend = iesm_ext(isma)
      if (ismb == isma) ibsta = ibsta+1
      nna = nn0+ibsta-1
      do ib=ibsta,ibend
        nna = nna+1
        nnb = nn0+iasta-1
        valuelptmp1 = vector1(nna)
        valuelptmp2 = vector2(nna)
        do ia=iasta,min(iaend,ib-1)
          valuelp1 = value_lpext(ii)
          valuelp2 = value_lpext(ii+1)
          mm = mm+1
          nnb = nnb+1
          vector2(mm) = vector2(mm)+valuelptmp1*valuelp1+vector1(nnb)*valuelp2
          vector2(nnb) = vector2(nnb)+vector1(mm)*valuelp2
          valuelptmp2 = valuelptmp2+vector1(mm)*valuelp1
          ii = ii+2

        end do
        vector2(nna) = valuelptmp2
      end do
    end do
  end if
  if (logic_g2g4b) then
    ii = ii0
    mm0 = iplprwei
    nn0 = iplplwei+iwt_sm_s_ext
    mm = mm0             !severe_error
    do ismb=1,ng_sm
      isma = Mul(ismb,ism_g2g4)
      if (isma > ismb) cycle
      ibsta = ibsm_ext(ismb)
      ibend = iesm_ext(ismb)
      iasta = ibsm_ext(isma)
      iaend = iesm_ext(isma)
      if (ismb == isma) ibsta = ibsta+1
      nna = nn0+ibsta-1
      do ib=ibsta,ibend
        nna = nna+1
        nnb = nn0+iasta-1
        valuelptmp1 = vector1(nna)
        valuelptmp2 = vector2(nna)
        do ia=iasta,min(iaend,ib-1)
          valuelp2 = value_lpext(ii)
          valuelp1 = value_lpext(ii+1)
          mm = mm+1
          nnb = nnb+1
          vector2(mm) = vector2(mm)+valuelptmp1*valuelp1+vector1(nnb)*valuelp2
          vector2(nnb) = vector2(nnb)+vector1(mm)*valuelp2
          valuelptmp2 = valuelptmp2+vector1(mm)*valuelp1
          ii = ii+2
        end do
        vector2(nna) = valuelptmp2
      end do
    end do
  end if
  ii0 = ii-1         !severe_error

  do icle=1,2
    if (((icle == 1) .and. logic_g36a) .or. ((icle == 2) .and. logic_g36b)) then
      if (icle == 1) then
        mm0 = iplplwei
        nn0 = iplprwei
        lpsta36 = lpsta36a
        lpend36 = lpend36a
      else
        mm0 = iplprwei
        nn0 = iplplwei
        lpsta36 = lpsta36b
        lpend36 = lpend36b
      end if
      do iii=lpsta36,lpend36,4
        ilwtmp = lpext_wei(iii)
        irwtmp = lpext_wei(iii+1)
        ii = ii0+lpext_wei(iii+2)
        mm = mm0+ilwtmp
        nn = nn0+irwtmp
        do i=1,lpext_wei(iii+3)
          vector2(mm) = vector2(mm)+vector1(nn)*value_lpext(ii)
          vector2(nn) = vector2(nn)+vector1(mm)*value_lpext(ii)
          mm = mm+1
          nn = nn+1
        end do
      end do
    end if

    if ((icle == 1) .and. logic_g35a) then
      mm0 = iplplwei
      nn0 = iplprwei
      lpsta35 = lpsta35a
      lpend35 = lpend35a
      do iii=lpsta35,lpend35,4
        mm = mm0+lpext_wei(iii)
        nn = nn0+lpext_wei(iii+1)
        iij = ii0+lpext_wei(iii+2)
        valuelp = vector2(nn)
        valuelptmp1 = vector1(nn)
        do i=1,lpext_wei(iii+3)
          valuetmp = -value_lpext(iij)
          vector2(mm) = vector2(mm)+valuelptmp1*valuetmp
          valuelp = valuelp+vector1(mm)*valuetmp
          iij = iij+1
          mm = mm+1
        end do
        vector2(nn) = valuelp
      end do
    else if ((icle == 2) .and. logic_g35b) then
      mm0 = iplprwei
      nn0 = iplplwei
      lpsta35 = lpsta35b
      lpend35 = lpend35b
      do iii=lpsta35,lpend35,4
        mm = mm0+lpext_wei(iii)
        nn = nn0+lpext_wei(iii+1)
        iij = ii0+lpext_wei(iii+2)
        valuelp = vector2(nn)         !severe_error
        valuelptmp1 = vector1(nn)
        do i=1,lpext_wei(iii+3)
          vector2(mm) = vector2(mm)+valuelptmp1*value_lpext(iij)
          valuelp = valuelp+vector1(mm)*value_lpext(iij)
          iij = iij+1
          mm = mm+1
        end do
        vector2(nn) = valuelp      !severe_error
      end do
    end if

    if (((icle == 1) .and. logic_g34a) .or. ((icle == 2) .and. logic_g34b)) then
      if (icle == 1) then
        mm0 = iplplwei
        nn0 = iplprwei
        lpsta34 = lpsta34a
        lpend34 = lpend34a
      else
        mm0 = iplprwei
        nn0 = iplplwei
        lpsta34 = lpsta34b
        lpend34 = lpend34b
      end if
      do iii=lpsta34,lpend34,4
        ilwtmp = lpext_wei(iii)
        irwtmp = lpext_wei(iii+1)
        mm = mm0+ilwtmp
        nn = nn0+irwtmp
        iij = ii0+lpext_wei(iii+2)
        valuelp = vector2(nn)
        valuelptmp1 = vector1(nn)
        do i=1,lpext_wei(iii+3)
          valuetmp = -value_lpext(iij)
          vector2(mm) = vector2(mm)+valuelptmp1*valuetmp
          valuelp = valuelp+vector1(mm)*valuetmp
          iij = iij+1
          mm = mm+1
        end do
        vector2(nn) = valuelp      !severe_error
      end do
    end if
    ii0 = ii0+nvalue_space_ss
  end do
end do

end subroutine inn_ext_ts_loop_unpack

subroutine inn_ext_ts_drl_loop_unpack(iplplweiin,iplprweiin)
!***********************************************************************
! 26 feb 2007 - revised

use gugaci_global, only: ibsm_ext, idownwei_g131415, iesm_ext, indx, ism_g2g4, iwt_sm_s_ext, logic_g1415, logic_g2g4a, &
                         logic_g2g4b, logic_g34a, logic_g34b, logic_g35a, logic_g35b, logic_g36a, logic_g36b, lpend34a, lpend34b, &
                         lpend35a, lpend35b, lpend36a, lpend36b, lpext_wei, lpsta34a, lpsta34b, lpsta35a, lpsta35b, lpsta36a, &
                         lpsta36b, mcroot, ng_sm, value_lpext, vector1, vector2
use Symmetry_Info, only: Mul
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: iplplweiin, iplprweiin
integer(kind=iwp) :: i, ia, iaend, iasta, ib, ibend, ibsta, icle, ii, ii0, iii, iij, ilwtmp, iplplwei, iplprwei, irot, irtidx, &
                     irwtmp, isma, ismb, lpend34, lpend35, lpend36, lpsta34, lpsta35, lpsta36, mm, mm0, nn, nn0, nna, nnb
real(kind=wp) :: valuelp, valuelp1, valuelp2, valuelptmp1, valuelptmp2, valuetmp

!write(u6,*) '  ts_test 2/2'

do irot=1,mcroot
  irtidx = indx(irot)
  iplplwei = iplplweiin+irtidx
  iplprwei = iplprweiin+irtidx

  ii = 1
  if (logic_g1415) then
    mm = iplplwei
    nn = iplprwei
    do i=1,idownwei_g131415
      valuelp = value_lpext(ii)
      mm = mm+1
      nn = nn+1
      vector2(mm) = vector2(mm)+vector1(nn)*valuelp
      vector2(nn) = vector2(nn)+vector1(mm)*valuelp
      ii = ii+1
    end do
  end if

  ii0 = ii
  if (logic_g2g4a) then
    ii = ii0
    mm0 = iplplwei
    nn0 = iplprwei+iwt_sm_s_ext
    mm = mm0            !severe_error
    do ismb=1,ng_sm
      isma = Mul(ismb,ism_g2g4)
      if (isma > ismb) cycle
      ibsta = ibsm_ext(ismb)
      ibend = iesm_ext(ismb)
      iasta = ibsm_ext(isma)
      iaend = iesm_ext(isma)
      if (ismb == isma) ibsta = ibsta+1
      nna = nn0+ibsta-1
      do ib=ibsta,ibend
        nna = nna+1
        nnb = nn0+iasta-1
        valuelptmp1 = vector1(nna)
        valuelptmp2 = vector2(nna)
        do ia=iasta,min(iaend,ib-1)
          valuelp1 = value_lpext(ii)
          valuelp2 = value_lpext(ii+1)
          mm = mm+1
          nnb = nnb+1
          vector2(mm) = vector2(mm)+valuelptmp1*valuelp1+vector1(nnb)*valuelp2
          vector2(nnb) = vector2(nnb)+vector1(mm)*valuelp2
          valuelptmp2 = valuelptmp2+vector1(mm)*valuelp1
          ii = ii+2
        end do
        vector2(nna) = valuelptmp2
      end do
    end do
  end if

  if (logic_g2g4b) then
    ii = ii0
    mm0 = iplprwei
    nn0 = iplplwei+iwt_sm_s_ext
    mm = mm0             !severe_error
    do ismb=1,ng_sm
      isma = Mul(ismb,ism_g2g4)
      if (isma > ismb) cycle
      ibsta = ibsm_ext(ismb)
      ibend = iesm_ext(ismb)
      iasta = ibsm_ext(isma)
      iaend = iesm_ext(isma)
      if (ismb == isma) ibsta = ibsta+1
      nna = nn0+ibsta-1
      do ib=ibsta,ibend
        nna = nna+1
        nnb = nn0+iasta-1
        valuelptmp1 = vector1(nna)
        valuelptmp2 = vector2(nna)
        do ia=iasta,min(iaend,ib-1)
          valuelp2 = value_lpext(ii)
          valuelp1 = value_lpext(ii+1)
          mm = mm+1
          nnb = nnb+1
          vector2(mm) = vector2(mm)+valuelptmp1*valuelp1+vector1(nnb)*valuelp2
          vector2(nnb) = vector2(nnb)+vector1(mm)*valuelp2
          valuelptmp2 = valuelptmp2+vector1(mm)*valuelp1
          ii = ii+2
        end do
        vector2(nna) = valuelptmp2
      end do
    end do
  end if

  ii0 = ii-1         !severe_error

  do icle=1,2
    if (((icle == 1) .and. logic_g36a) .or. ((icle == 2) .and. logic_g36b)) then
      if (icle == 1) then
        mm0 = iplplwei
        nn0 = iplprwei
        lpsta36 = lpsta36a
        lpend36 = lpend36a
      else
        mm0 = iplprwei
        nn0 = iplplwei
        lpsta36 = lpsta36b
        lpend36 = lpend36b
      end if
      do iii=lpsta36,lpend36,4
        ilwtmp = lpext_wei(iii)
        irwtmp = lpext_wei(iii+1)
        ii = ii0+lpext_wei(iii+2)
        mm = mm0+ilwtmp
        nn = nn0+irwtmp
        do i=1,lpext_wei(iii+3)
          vector2(mm) = vector2(mm)+vector1(nn)*value_lpext(ii)
          vector2(nn) = vector2(nn)+vector1(mm)*value_lpext(ii)
          mm = mm+1
          nn = nn+1
        end do
      end do
    end if

    if ((icle == 1) .and. logic_g35a) then
      mm0 = iplplwei
      nn0 = iplprwei
      lpsta35 = lpsta35a
      lpend35 = lpend35a
      do iii=lpsta35,lpend35,4
        mm = mm0+lpext_wei(iii)
        nn = nn0+lpext_wei(iii+1)
        iij = ii0+lpext_wei(iii+2)
        valuelp = vector2(nn)
        valuelptmp1 = vector1(nn)
        do i=1,lpext_wei(iii+3)
          valuetmp = -value_lpext(iij)
          vector2(mm) = vector2(mm)+valuelptmp1*valuetmp
          valuelp = valuelp+vector1(mm)*valuetmp
          iij = iij+1
          mm = mm+1
        end do
        vector2(nn) = valuelp
      end do
    else if ((icle == 2) .and. logic_g35b) then
      mm0 = iplprwei
      nn0 = iplplwei
      lpsta35 = lpsta35b
      lpend35 = lpend35b
      do iii=lpsta35,lpend35,4
        mm = mm0+lpext_wei(iii)
        nn = nn0+lpext_wei(iii+1)
        iij = ii0+lpext_wei(iii+2)
        valuelp = vector2(nn)         !severe_error
        valuelptmp1 = vector1(nn)
        do i=1,lpext_wei(iii+3)
          vector2(mm) = vector2(mm)+valuelptmp1*value_lpext(iij)
          valuelp = valuelp+vector1(mm)*value_lpext(iij)
          iij = iij+1
          mm = mm+1
        end do
        vector2(nn) = valuelp      !severe_error
      end do
    end if

    if (((icle == 1) .and. logic_g34a) .or. ((icle == 2) .and. logic_g34b)) then
      if (icle == 1) then
        mm0 = iplplwei
        nn0 = iplprwei
        lpsta34 = lpsta34a
        lpend34 = lpend34a
      else
        mm0 = iplprwei
        nn0 = iplplwei
        lpsta34 = lpsta34b
        lpend34 = lpend34b
      end if
      do iii=lpsta34,lpend34,4
        ilwtmp = lpext_wei(iii)
        irwtmp = lpext_wei(iii+1)
        mm = mm0+ilwtmp
        nn = nn0+irwtmp
        iij = ii0+lpext_wei(iii+2)
        valuelp = vector2(nn)
        valuelptmp1 = vector1(nn)
        do i=1,lpext_wei(iii+3)
          valuetmp = -value_lpext(iij)
          vector2(mm) = vector2(mm)+valuelptmp1*valuetmp
          valuelp = valuelp+vector1(mm)*valuetmp
          iij = iij+1
          mm = mm+1
        end do
        vector2(nn) = valuelp      !severe_error
      end do
    end if
  end do

end do

end subroutine inn_ext_ts_drl_loop_unpack

subroutine inn_ext_st_loop_unpack(iplplweiin,iplprweiin)
!***********************************************************************
! 26 feb 2007 - revised

use gugaci_global, only: ibsm_ext, idownwei_g131415, iesm_ext, indx, ism_g2g4, iwt_sm_s_ext, logic_g1415, logic_g2g4a, &
                         logic_g2g4b, logic_g34a, logic_g34b, logic_g35a, logic_g35b, logic_g36a, logic_g36b, lpend34a, lpend34b, &
                         lpend35a, lpend35b, lpend36a, lpend36b, lpext_wei, lpsta34a, lpsta34b, lpsta35a, lpsta35b, lpsta36a, &
                         lpsta36b, mcroot, ng_sm, nvalue_space_ss, value_lpext, vector1, vector2
use Symmetry_Info, only: Mul
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: iplplweiin, iplprweiin
integer(kind=iwp) :: i, ia, iaend, iasta, ib, ibend, ibsta, icle, ii, ii0, iii, iij, ilwtmp, iplplwei, iplprwei, irot, irtidx, &
                     irwtmp, isma, ismb, lpend34, lpend35, lpend36, lpsta34, lpsta35, lpsta36, mm, mm0, nn, nn0, nna, nnb
real(kind=wp) :: valtmp, valuelp, valuelp1, valuelp2, valuelptmp1, valuelptmp2, valuetmp

!write(u6,*) ' st_test 1/2'

do irot=1,mcroot
  irtidx = indx(irot)
  iplplwei = iplplweiin+irtidx
  iplprwei = iplprweiin+irtidx

  ii = 1
  if (logic_g1415) then
    mm = iplplwei
    nn = iplprwei
    do i=1,idownwei_g131415
      valuelp = value_lpext(ii)
      mm = mm+1
      nn = nn+1
      vector2(mm) = vector2(mm)+vector1(nn)*valuelp
      vector2(nn) = vector2(nn)+vector1(mm)*valuelp
      ii = ii+1
    end do
  end if

  !write(u6,*) 'st_g2g4a',iplplwei,iplprwei,vector2(137)
  ii0 = ii
  if (logic_g2g4a) then
    ii = ii0
    mm0 = iplplwei
    nn0 = iplprwei+iwt_sm_s_ext
    mm = mm0            !severe_error
    do ismb=1,ng_sm
      isma = Mul(ismb,ism_g2g4)
      if (isma > ismb) cycle
      ibsta = ibsm_ext(ismb)
      ibend = iesm_ext(ismb)
      iasta = ibsm_ext(isma)
      iaend = iesm_ext(isma)
      if (ismb == isma) ibsta = ibsta+1
      nna = nn0+ibsta-1
      do ib=ibsta,ibend
        nna = nna+1
        nnb = nn0+iasta-1
        valuelptmp1 = vector1(nna)
        valuelptmp2 = vector2(nna)
        do ia=iasta,min(iaend,ib-1)
          valuelp1 = value_lpext(ii)
          valuelp2 = value_lpext(ii+1)
          mm = mm+1
          nnb = nnb+1
          vector2(mm) = vector2(mm)+valuelptmp1*valuelp1+vector1(nnb)*valuelp2
          vector2(nnb) = vector2(nnb)+vector1(mm)*valuelp2
          valuelptmp2 = valuelptmp2+vector1(mm)*valuelp1
          ii = ii+2
        end do
        vector2(nna) = valuelptmp2
      end do
    end do
  end if
  !write(u6,*) 'st_g2g4b',iplplwei,iplprwei,vector2(137)
  if (logic_g2g4b) then
    ii = ii0
    mm0 = iplprwei
    nn0 = iplplwei+iwt_sm_s_ext
    mm = mm0             !severe_error
    do ismb=1,ng_sm
      isma = Mul(ismb,ism_g2g4)
      if (isma > ismb) cycle
      ibsta = ibsm_ext(ismb)
      ibend = iesm_ext(ismb)
      iasta = ibsm_ext(isma)
      iaend = iesm_ext(isma)
      if (ismb == isma) ibsta = ibsta+1
      nna = nn0+ibsta-1
      do ib=ibsta,ibend
        nna = nna+1
        nnb = nn0+iasta-1
        valuelptmp1 = vector1(nna)
        valuelptmp2 = vector2(nna)
        do ia=iasta,min(iaend,ib-1)
          valuelp2 = value_lpext(ii)
          valuelp1 = value_lpext(ii+1)
          mm = mm+1
          nnb = nnb+1
          vector2(mm) = vector2(mm)+valuelptmp1*valuelp1+vector1(nnb)*valuelp2
          vector2(nnb) = vector2(nnb)+vector1(mm)*valuelp2
          valuelptmp2 = valuelptmp2+vector1(mm)*valuelp1
          ii = ii+2
        end do
        vector2(nna) = valuelptmp2
      end do
    end do
  end if
  ii0 = ii-1

  !write(u6,*) 'st_g36a',iplplwei,iplprwei,vector2(137)
  do icle=1,2
    if (((icle == 1) .and. logic_g36a) .or. ((icle == 2) .and. logic_g36b)) then
      if (icle == 1) then
        mm0 = iplplwei
        nn0 = iplprwei
        lpsta36 = lpsta36a
        lpend36 = lpend36a
      else
        mm0 = iplprwei
        nn0 = iplplwei
        lpsta36 = lpsta36b
        lpend36 = lpend36b
      end if
      do iii=lpsta36,lpend36,4
        ilwtmp = lpext_wei(iii)
        irwtmp = lpext_wei(iii+1)
        ii = ii0+lpext_wei(iii+2)
        mm = mm0+ilwtmp
        nn = nn0+irwtmp
        do i=1,lpext_wei(iii+3)
          vector2(mm) = vector2(mm)+vector1(nn)*value_lpext(ii)
          vector2(nn) = vector2(nn)+vector1(mm)*value_lpext(ii)
          mm = mm+1
          nn = nn+1
        end do
      end do
    end if

    !write(u6,*) 'st_g35a',iplplwei,iplprwei,vector2(137)
    if ((icle == 1) .and. logic_g35a) then
      mm0 = iplplwei
      nn0 = iplprwei
      lpsta35 = lpsta35a
      lpend35 = lpend35a
      do iii=lpsta35,lpend35,4
        mm = mm0+lpext_wei(iii)
        nn = nn0+lpext_wei(iii+1)
        iij = ii0+lpext_wei(iii+2)
        valuelp = vector2(nn)
        valuelptmp1 = vector1(nn)
        do i=1,lpext_wei(iii+3)
          valtmp = value_lpext(iij)
          vector2(mm) = vector2(mm)+valuelptmp1*valtmp
          valuelp = valuelp+vector1(mm)*valtmp
          iij = iij+1
          mm = mm+1
        end do
        vector2(nn) = valuelp
      end do
    else if ((icle == 2) .and. logic_g35b) then
      mm0 = iplprwei
      nn0 = iplplwei
      lpsta35 = lpsta35b
      lpend35 = lpend35b
      do iii=lpsta35,lpend35,4
        mm = mm0+lpext_wei(iii)
        nn = nn0+lpext_wei(iii+1)
        iij = ii0+lpext_wei(iii+2)
        valuelp = vector2(nn)
        valuelptmp1 = vector1(nn)
        do i=1,lpext_wei(iii+3)
          valtmp = -value_lpext(iij)
          vector2(mm) = vector2(mm)+valuelptmp1*valtmp
          valuelp = valuelp+vector1(mm)*valtmp
          iij = iij+1
          mm = mm+1
        end do
        vector2(nn) = valuelp      !severe_error
      end do
    end if

    !write(u6,*) 'st_g34',iplplwei,iplprwei,vector2(137)
    if (((icle == 1) .and. logic_g34a) .or. ((icle == 2) .and. logic_g34b)) then
      if (icle == 1) then
        mm0 = iplplwei
        nn0 = iplprwei
        lpsta34 = lpsta34a
        lpend34 = lpend34a
      else
        mm0 = iplprwei
        nn0 = iplplwei
        lpsta34 = lpsta34b
        lpend34 = lpend34b
      end if
      do iii=lpsta34,lpend34,4
        ilwtmp = lpext_wei(iii)
        irwtmp = lpext_wei(iii+1)
        mm = mm0+ilwtmp
        nn = nn0+irwtmp
        iij = ii0+lpext_wei(iii+2)
        valuelp = vector2(nn)
        valuelptmp1 = vector1(nn)
        do i=1,lpext_wei(iii+3)
          valuetmp = -value_lpext(iij)
          vector2(mm) = vector2(mm)+valuelptmp1*valuetmp
          valuelp = valuelp+vector1(mm)*valuetmp
          iij = iij+1
          mm = mm+1
        end do
        vector2(nn) = valuelp
      end do
    end if
    ii0 = ii0+nvalue_space_ss
  end do

end do

end subroutine inn_ext_st_loop_unpack

subroutine inn_ext_st_drl_loop_unpack(iplplweiin,iplprweiin)
!***********************************************************************
! 26 feb 2007 - revised

use gugaci_global, only: ibsm_ext, idownwei_g131415, iesm_ext, indx, ism_g2g4, iwt_sm_s_ext, logic_g1415, logic_g2g4a, &
                         logic_g2g4b, logic_g34a, logic_g34b, logic_g35a, logic_g35b, logic_g36a, logic_g36b, lpend34a, lpend34b, &
                         lpend35a, lpend35b, lpend36a, lpend36b, lpext_wei, lpsta34a, lpsta34b, lpsta35a, lpsta35b, lpsta36a, &
                         lpsta36b, mcroot, ng_sm, value_lpext, vector1, vector2
use Symmetry_Info, only: Mul
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: iplplweiin, iplprweiin
integer(kind=iwp) :: i, ia, iaend, iasta, ib, ibend, ibsta, icle, ii, ii0, iii, iij, ilwtmp, iplplwei, iplprwei, irot, irtidx, &
                     irwtmp, isma, ismb, lpend34, lpend35, lpend36, lpsta34, lpsta35, lpsta36, mm, mm0, nn, nn0, nna, nnb
real(kind=wp) :: valuelp, valuelp1, valuelp2, valuelptmp1, valuelptmp2, valuetmp

!write(u6,*) ' st_test 2/2'

do irot=1,mcroot
  irtidx = indx(irot)
  iplplwei = iplplweiin+irtidx
  iplprwei = iplprweiin+irtidx

  ii = 1
  if (logic_g1415) then
    mm = iplplwei
    nn = iplprwei
    do i=1,idownwei_g131415
      valuelp = value_lpext(ii)
      mm = mm+1
      nn = nn+1
      vector2(mm) = vector2(mm)+vector1(nn)*valuelp
      vector2(nn) = vector2(nn)+vector1(mm)*valuelp
      ii = ii+1
    end do
  end if

  ii0 = ii
  if (logic_g2g4a) then
    ii = ii0
    mm0 = iplplwei
    nn0 = iplprwei+iwt_sm_s_ext
    mm = mm0            !severe_error
    do ismb=1,ng_sm
      isma = Mul(ismb,ism_g2g4)
      if (isma > ismb) cycle
      ibsta = ibsm_ext(ismb)
      ibend = iesm_ext(ismb)
      iasta = ibsm_ext(isma)
      iaend = iesm_ext(isma)
      if (ismb == isma) ibsta = ibsta+1
      nna = nn0+ibsta-1
      do ib=ibsta,ibend
        nna = nna+1
        nnb = nn0+iasta-1
        valuelptmp1 = vector1(nna)
        valuelptmp2 = vector2(nna)
        do ia=iasta,min(iaend,ib-1)
          valuelp1 = value_lpext(ii)
          valuelp2 = value_lpext(ii+1)
          mm = mm+1
          nnb = nnb+1
          vector2(mm) = vector2(mm)+valuelptmp1*valuelp1+vector1(nnb)*valuelp2
          vector2(nnb) = vector2(nnb)+vector1(mm)*valuelp2
          valuelptmp2 = valuelptmp2+vector1(mm)*valuelp1
          ii = ii+2
        end do
        vector2(nna) = valuelptmp2
      end do
    end do
  end if
  if (logic_g2g4b) then
    ii = ii0
    mm0 = iplprwei
    nn0 = iplplwei+iwt_sm_s_ext
    mm = mm0             !severe_error
    do ismb=1,ng_sm
      isma = Mul(ismb,ism_g2g4)
      if (isma > ismb) cycle
      ibsta = ibsm_ext(ismb)
      ibend = iesm_ext(ismb)
      iasta = ibsm_ext(isma)
      iaend = iesm_ext(isma)
      if (ismb == isma) ibsta = ibsta+1
      nna = nn0+ibsta-1
      do ib=ibsta,ibend
        nna = nna+1
        nnb = nn0+iasta-1
        valuelptmp1 = vector1(nna)
        valuelptmp2 = vector2(nna)
        do ia=iasta,min(iaend,ib-1)
          valuelp2 = value_lpext(ii)
          valuelp1 = value_lpext(ii+1)
          mm = mm+1
          nnb = nnb+1
          vector2(mm) = vector2(mm)+valuelptmp1*valuelp1+vector1(nnb)*valuelp2
          vector2(nnb) = vector2(nnb)+vector1(mm)*valuelp2
          valuelptmp2 = valuelptmp2+vector1(mm)*valuelp1
          ii = ii+2
        end do
        vector2(nna) = valuelptmp2
      end do
    end do
  end if
  ii0 = ii-1

  do icle=1,2
    if (((icle == 1) .and. logic_g36a) .or. ((icle == 2) .and. logic_g36b)) then
      if (icle == 1) then
        mm0 = iplplwei
        nn0 = iplprwei
        lpsta36 = lpsta36a
        lpend36 = lpend36a
      else
        mm0 = iplprwei
        nn0 = iplplwei
        lpsta36 = lpsta36b
        lpend36 = lpend36b
      end if
      do iii=lpsta36,lpend36,4
        ilwtmp = lpext_wei(iii)
        irwtmp = lpext_wei(iii+1)
        ii = ii0+lpext_wei(iii+2)
        mm = mm0+ilwtmp
        nn = nn0+irwtmp
        do i=1,lpext_wei(iii+3)
          valuetmp = value_lpext(ii)
          vector2(mm) = vector2(mm)+vector1(nn)*valuetmp
          vector2(nn) = vector2(nn)+vector1(mm)*valuetmp
          mm = mm+1
          nn = nn+1
        end do
      end do
    end if
    if ((icle == 1) .and. logic_g35a) then
      mm0 = iplplwei
      nn0 = iplprwei
      lpsta35 = lpsta35a
      lpend35 = lpend35a
      do iii=lpsta35,lpend35,4
        mm = mm0+lpext_wei(iii)
        nn = nn0+lpext_wei(iii+1)
        iij = ii0+lpext_wei(iii+2)
        valuelp = vector2(nn)
        valuelptmp1 = vector1(nn)
        do i=1,lpext_wei(iii+3)
          valuetmp = value_lpext(iij)
          vector2(mm) = vector2(mm)+valuelptmp1*valuetmp
          valuelp = valuelp+vector1(mm)*valuetmp
          iij = iij+1
          mm = mm+1
        end do
        vector2(nn) = valuelp
      end do
    else if ((icle == 2) .and. logic_g35b) then
      mm0 = iplprwei
      nn0 = iplplwei
      lpsta35 = lpsta35b
      lpend35 = lpend35b
      do iii=lpsta35,lpend35,4
        mm = mm0+lpext_wei(iii)
        nn = nn0+lpext_wei(iii+1)
        iij = ii0+lpext_wei(iii+2)
        valuelp = vector2(nn)         !severe_error
        valuelptmp1 = vector1(nn)
        do i=1,lpext_wei(iii+3)
          valuetmp = -value_lpext(iij)
          vector2(mm) = vector2(mm)+valuelptmp1*valuetmp
          valuelp = valuelp+vector1(mm)*valuetmp
          iij = iij+1
          mm = mm+1
        end do
        vector2(nn) = valuelp      !severe_error
      end do
    end if

    if (((icle == 1) .and. logic_g34a) .or. ((icle == 2) .and. logic_g34b)) then
      if (icle == 1) then
        mm0 = iplplwei
        nn0 = iplprwei
        lpsta34 = lpsta34a
        lpend34 = lpend34a
      else
        mm0 = iplprwei
        nn0 = iplplwei
        lpsta34 = lpsta34b
        lpend34 = lpend34b
      end if
      do iii=lpsta34,lpend34,4
        ilwtmp = lpext_wei(iii)
        irwtmp = lpext_wei(iii+1)
        mm = mm0+ilwtmp
        nn = nn0+irwtmp
        iij = ii0+lpext_wei(iii+2)
        valuelp = vector2(nn)
        valuelptmp1 = vector1(nn)
        do i=1,lpext_wei(iii+3)
          valuetmp = -value_lpext(iij)
          vector2(mm) = vector2(mm)+valuelptmp1*valuetmp
          valuelp = valuelp+vector1(mm)*valuetmp
          iij = iij+1
          mm = mm+1
        end do
        vector2(nn) = valuelp
      end do
    end if
  end do
end do

end subroutine inn_ext_st_drl_loop_unpack

subroutine gsd_sequence_extspace(iplplweiin,iplprweiin)
!***********************************************************************
! 26 feb 2007 - revised

use gugaci_global, only: indx, ivaluesta_g26, iweista_g25, iweista_g26, iweista_g28, logic_g25a, logic_g25b, logic_g26, &
                         logic_g28a, mcroot, nint_g25, nint_g28, nwei_g25, nwei_g26, nwei_g28, v_sqtwo, value_lpext, vector1, &
                         vector2
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: iplplweiin, iplprweiin
integer(kind=iwp) :: i, ilpvalue, iplplwei, iplprwei, irot, irtidx, itmp, mm, nn, nn0
real(kind=wp) :: vlptmp, vlptmp1, vtmp

!write(u6,*) '  sd_test 1/2','  ds_test 011'

do irot=1,mcroot
  irtidx = indx(irot)
  iplplwei = iplplweiin+irtidx
  iplprwei = iplprweiin+irtidx
  ilpvalue = 0
  if (logic_g25a) then
    mm = iplplwei+iweista_g25-1
    nn0 = iplprwei
    do itmp=1,nint_g25
      ilpvalue = ilpvalue+1
      vlptmp = value_lpext(ilpvalue)
      nn = nn0
      do i=1,nwei_g25
        mm = mm+1
        nn = nn+1
        vector2(mm) = vector2(mm)+vector1(nn)*vlptmp
        vector2(nn) = vector2(nn)+vector1(mm)*vlptmp
      end do
    end do
  else if (logic_g25b) then
    mm = iplplwei+iweista_g25-1
    nn0 = iplprwei
    ilpvalue = ilpvalue+1
    do itmp=2,nint_g25
      ilpvalue = ilpvalue+1
      vlptmp = value_lpext(ilpvalue)
      nn = nn0
      do i=1,itmp-1
        mm = mm+1
        nn = nn+1
        vector2(mm) = vector2(mm)+vector1(nn)*vlptmp
        vector2(nn) = vector2(nn)+vector1(mm)*vlptmp
      end do
    end do

    mm = iplplwei+iweista_g28-1
    nn = iplprwei
    nn = nn+1
    do itmp=2,nwei_g28
      nn = nn+1
      vlptmp = vector2(nn)
      vlptmp1 = vector1(nn)
      ilpvalue = 0
      do i=1,itmp-1
        ilpvalue = ilpvalue+1
        vtmp = value_lpext(ilpvalue)
        mm = mm+1
        vector2(mm) = vector2(mm)+vlptmp1*vtmp
        vlptmp = vlptmp+vector1(mm)*vtmp
      end do
      vector2(nn) = vlptmp
    end do
  else if (logic_g28a) then
    mm = iplplwei+iweista_g28-1
    nn0 = iplprwei
    do nn=nn0+1,nn0+nwei_g28
      vlptmp = vector2(nn)
      vlptmp1 = vector1(nn)
      ilpvalue = 0
      do i=1,nint_g28
        ilpvalue = ilpvalue+1
        mm = mm+1
        vector2(mm) = vector2(mm)+vlptmp1*value_lpext(ilpvalue)
        vlptmp = vlptmp+vector1(mm)*value_lpext(ilpvalue)
      end do
      vector2(nn) = vlptmp
    end do
  end if

  if (logic_g26) then
    ilpvalue = ivaluesta_g26
    mm = iplplwei+iweista_g26
    nn0 = iplprwei
    do nn=nn0+1,nn0+nwei_g26
      ilpvalue = ilpvalue+1
      vlptmp = value_lpext(ilpvalue)*v_sqtwo   !for g26
      vector2(mm) = vector2(mm)+vector1(nn)*vlptmp
      vector2(nn) = vector2(nn)+vector1(mm)*vlptmp
      mm = mm+1
    end do
  end if
end do

end subroutine gsd_sequence_extspace

subroutine complete_sd_ar_ext_loop(ilweiin,irweiin,isdownwei)
!***********************************************************************
! 26 feb 2007 - revised

use gugaci_global, only: icano_nnend, icano_nnsta, indx, mcroot, value_lpext, vector1, vector2
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: ilweiin, irweiin, isdownwei
integer(kind=iwp) :: ilpvalue, ilwei, irot, irtidx, irwei, mm, mm0, mmtmp, nn, nntmp
real(kind=wp) :: vlptmp, vlptmp1

!write(u6,*) 'sd_test 2/2','  td_test_2/2 012'

do irot=1,mcroot
  irtidx = indx(irot)
  ilwei = ilweiin+irtidx
  irwei = irweiin+irtidx

  ilpvalue = 0
  mm0 = ilwei
  nn = irwei+icano_nnsta-1
  do nntmp=icano_nnsta,icano_nnend
    nn = nn+1
    mm = mm0
    vlptmp1 = vector1(nn)
    vlptmp = vector2(nn)
    do mmtmp=1,isdownwei
      ilpvalue = ilpvalue+1
      mm = mm+1
      vector2(mm) = vector2(mm)+vlptmp1*value_lpext(ilpvalue)
      vlptmp = vlptmp+vector1(mm)*value_lpext(ilpvalue)
    end do
    vector2(nn) = vlptmp
  end do
end do

return

end subroutine complete_sd_ar_ext_loop

subroutine inn_ext_dd_loop_unpack(iplplweiin,iplprweiin)
!***********************************************************************
! 26 feb 2007 - revised

use gugaci_global, only: ildownwei_segdd, indx, int_dd_drl, irdownwei_segdd, logic_g49a, logic_g49b, logic_g50, mcroot, &
                         value_lpext, vector1, vector2
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: iplplweiin, iplprweiin
integer(kind=iwp) :: i, icle, ii, ildownwei, iplplwei, iplprwei, irdownwei, irot, irtidx, j, mm, mm0, nn
real(kind=wp) :: valuelp, valuelptmp1, valuetmp

! 'dd_test '

do irot=1,mcroot
  irtidx = indx(irot)
  iplplwei = iplplweiin+irtidx
  iplprwei = iplprweiin+irtidx

  ii = 1
  if (logic_g50) then         !arbl=.true.
    if (logic_g49b) then
      mm = iplplwei
      nn = iplprwei
      do i=1,ildownwei_segdd
        valuelp = value_lpext(ii)
        mm = mm+1
        nn = nn+1
        vector2(mm) = vector2(mm)+vector1(nn)*valuelp
        vector2(nn) = vector2(nn)+vector1(mm)*valuelp
        ii = ii+1
      end do
    end if

    ii = ii+int_dd_drl
    mm0 = iplplwei
    nn = iplprwei+1
    do icle=1,2
      do j=2,ildownwei_segdd
        nn = nn+1
        valuelp = vector2(nn)
        valuelptmp1 = vector1(nn)
        mm = mm0
        do i=1,j-1
          mm = mm+1
          valuetmp = value_lpext(ii)
          vector2(mm) = vector2(mm)+valuelptmp1*valuetmp
          valuelp = valuelp+vector1(mm)*valuetmp
          ii = ii+1
        end do
        vector2(nn) = valuelp
      end do
      if (.not. logic_g49b) exit
      mm0 = iplprwei
      nn = iplplwei+1
    end do

  else               !drl=.true.

    ii = ii+int_dd_drl
    if (logic_g49a) then
      mm0 = iplplwei
      nn = iplprwei
      ildownwei = ildownwei_segdd
      irdownwei = irdownwei_segdd
    else
      mm0 = iplprwei
      nn = iplplwei
      ildownwei = irdownwei_segdd
      irdownwei = ildownwei_segdd
    end if
    do j=1,irdownwei
      nn = nn+1
      valuelp = vector2(nn)
      valuelptmp1 = vector1(nn)
      mm = mm0
      do i=1,ildownwei
        mm = mm+1
        valuetmp = value_lpext(ii)
        vector2(mm) = vector2(mm)+valuelptmp1*valuetmp
        valuelp = valuelp+vector1(mm)*valuetmp
        ii = ii+1
      end do
      vector2(nn) = valuelp
    end do
  end if
end do

end subroutine inn_ext_dd_loop_unpack

subroutine inn_ext_tt_drl_loop_unpack(iplplweiin,iplprweiin,n1415)
!***********************************************************************
! 26 feb 2007 - revised

use gugaci_global, only: idownwei_g131415, indx, logic_g1415, logic_g34a, logic_g34b, logic_g35a, logic_g35b, logic_g36a, &
                         logic_g36b, lpend34a, lpend34b, lpend35a, lpend35b, lpend36a, lpend36b, lpext_wei, lpsta34a, lpsta34b, &
                         lpsta35a, lpsta35b, lpsta36a, lpsta36b, mcroot, value_lpext, vector1, vector2
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: iplplweiin, iplprweiin, n1415
integer(kind=iwp) :: i, icle, ii, ii0, iii, iij, ilwtmp, iplplwei, iplprwei, irot, irtidx, irwtmp, lpend34, lpend35, lpend36, &
                     lpsta34, lpsta35, lpsta36, mm, mm0, nn, nn0
real(kind=wp) :: valuelp, valuelptmp1, valuetmp
logical(kind=iwp) :: logic_g14150, logic_g34b0, logic_g35b0, logic_g36b0

!write(u6,*) '  tt_test 2/2'

logic_g14150 = logic_g1415
logic_g36b0 = logic_g36b
logic_g35b0 = logic_g35b
logic_g34b0 = logic_g34b
if (iplplweiin == iplprweiin) then
  logic_g14150 = .false.
  logic_g36b0 = .false.
  logic_g35b0 = .false.
  logic_g34b0 = .false.
end if

do irot=1,mcroot
  irtidx = indx(irot)
  iplplwei = iplplweiin+irtidx
  iplprwei = iplprweiin+irtidx

  ii = 1
  if (logic_g14150) then
    mm = iplplwei
    nn = iplprwei
    do i=1,idownwei_g131415
      valuelp = value_lpext(ii)
      mm = mm+1
      nn = nn+1
      vector2(mm) = vector2(mm)+vector1(nn)*valuelp
      vector2(nn) = vector2(nn)+vector1(mm)*valuelp
      ii = ii+1
    end do
  end if

  ii0 = n1415         !severe_error

  do icle=1,2
    if (((icle == 1) .and. logic_g36a) .or. ((icle == 2) .and. logic_g36b0)) then
      if (icle == 1) then
        mm0 = iplplwei
        nn0 = iplprwei
        lpsta36 = lpsta36a
        lpend36 = lpend36a
      else
        mm0 = iplprwei
        nn0 = iplplwei
        lpsta36 = lpsta36b
        lpend36 = lpend36b
      end if
      do iii=lpsta36,lpend36,4
        ilwtmp = lpext_wei(iii)
        irwtmp = lpext_wei(iii+1)
        ii = ii0+lpext_wei(iii+2)
        mm = mm0+ilwtmp
        nn = nn0+irwtmp
        do i=1,lpext_wei(iii+3)
          valuelp = value_lpext(ii)
          vector2(mm) = vector2(mm)+vector1(nn)*valuelp
          vector2(nn) = vector2(nn)+vector1(mm)*valuelp
          mm = mm+1
          nn = nn+1
        end do
      end do
    end if

    if (((icle == 1) .and. logic_g35a) .or. ((icle == 2) .and. logic_g35b0)) then
      if (icle == 1) then
        mm0 = iplplwei
        nn0 = iplprwei
        lpsta35 = lpsta35a
        lpend35 = lpend35a
      else
        mm0 = iplprwei
        nn0 = iplplwei
        lpsta35 = lpsta35b
        lpend35 = lpend35b
      end if
      do iii=lpsta35,lpend35,4
        mm = mm0+lpext_wei(iii)
        nn = nn0+lpext_wei(iii+1)
        iij = ii0+lpext_wei(iii+2)
        valuelp = vector2(nn)
        valuelptmp1 = vector1(nn)
        do i=1,lpext_wei(iii+3)
          valuetmp = -value_lpext(iij)
          vector2(mm) = vector2(mm)+valuelptmp1*valuetmp
          valuelp = valuelp+vector1(mm)*valuetmp
          iij = iij+1
          mm = mm+1
        end do
        vector2(nn) = valuelp      !severe_error_1202
      end do
    end if

    if (((icle == 1) .and. logic_g34a) .or. ((icle == 2) .and. logic_g34b0)) then
      if (icle == 1) then
        mm0 = iplplwei
        nn0 = iplprwei
        lpsta34 = lpsta34a
        lpend34 = lpend34a
      else
        mm0 = iplprwei
        nn0 = iplplwei
        lpsta34 = lpsta34b
        lpend34 = lpend34b
      end if
      do iii=lpsta34,lpend34,4
        ilwtmp = lpext_wei(iii)
        irwtmp = lpext_wei(iii+1)
        mm = mm0+ilwtmp
        nn = nn0+irwtmp
        iij = ii0+lpext_wei(iii+2)
        valuelp = vector2(nn)
        valuelptmp1 = vector1(nn)
        do i=1,lpext_wei(iii+3)
          vector2(mm) = vector2(mm)+valuelptmp1*value_lpext(iij)
          valuelp = valuelp+vector1(mm)*value_lpext(iij)
          iij = iij+1
          mm = mm+1
        end do
        vector2(nn) = valuelp
      end do
    end if
  end do
end do

end subroutine inn_ext_tt_drl_loop_unpack
