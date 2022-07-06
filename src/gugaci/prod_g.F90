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

subroutine inn_ext_ss_drl_loop_unpack_g(iplplwei,iplprwei)

use gugaci_global, only: ibsm_ext, iesm_ext, index_lpext, index_lpext1, ism_g2g4, iwt_sm_s_ext, logic_g2g4a, logic_g2g4b, &
                         logic_g34a, logic_g34b, logic_g35a, logic_g35b, logic_g36a, logic_g36b, lpend34a, lpend34b, lpend35a, &
                         lpend35b, lpend36a, lpend36b, lpext_wei, lpsta34a, lpsta34b, lpsta35a, lpsta35b, lpsta36a, lpsta36b, &
                         ng_sm, value_lpext, value_lpext1, vector1, vector2
use Symmetry_Info, only: Mul
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: iplplwei, iplprwei
integer(kind=iwp) :: i, ia, iaend, iasta, ib, ibend, ibsta, icle, ii, ii0, iii, iij, ilwtmp, indexlp, indexlp1, irwtmp, isma, &
                     ismb, lpend34, lpend35, lpend36, lpsta34, lpsta35, lpsta36, mm, mm0, nn, nn0, nna, nnb
real(kind=wp) :: valuelp, valuelp1

!write(u6,*) ' ss_test 2/2'
ii = 1
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

      do ia=iasta,min(iaend,ib-1)
        mm = mm+1
        nnb = nnb+1
        indexlp = index_lpext(ii+1)
        indexlp1 = index_lpext1(ii+1)

        !if (indexlp /= 0) then
        valuelp = value_lpext(ii+1)
        vector2(indexlp) = vector2(indexlp)+vector1(mm)*(vector1(nna)+vector1(nnb))*valuelp
        !end if
        !if (indexlp1 /= 0) then
        valuelp1 = value_lpext1(ii+1)
        vector2(indexlp1) = vector2(indexlp1)+vector1(mm)*(vector1(nna)+vector1(nnb))*valuelp1
        !end if
        ii = ii+2
      end do
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
      do ia=iasta,min(iaend,ib-1)
        mm = mm+1
        nnb = nnb+1
        indexlp = index_lpext(ii+1)
        indexlp1 = index_lpext1(ii+1)
        !if (indexlp /= 0) then
        valuelp = value_lpext(ii+1)
        vector2(indexlp) = vector2(indexlp)+vector1(mm)*(vector1(nna)+vector1(nnb))*valuelp
        !end if
        !if (indexlp1 /= 0) then
        valuelp1 = value_lpext1(ii+1)
        vector2(indexlp1) = vector2(indexlp1)+vector1(mm)*(vector1(nna)+vector1(nnb))*valuelp1
        !end if
        ii = ii+2
      end do
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
        indexlp = index_lpext(ii)
        indexlp1 = index_lpext1(ii)
        !if (indexlp /= 0) then
        valuelp = value_lpext(ii)
        vector2(indexlp) = vector2(indexlp)+vector1(mm)*vector1(nn)*valuelp
        !end if
        !if (indexlp1 /= 0) then
        valuelp1 = value_lpext1(ii)
        vector2(indexlp1) = vector2(indexlp1)+vector1(mm)*vector1(nn)*valuelp1
        !end if

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
      do i=1,lpext_wei(iii+3)
        indexlp = index_lpext(iij)
        indexlp1 = index_lpext1(iij)
        !if (indexlp /= 0) then
        valuelp = value_lpext(iij)
        vector2(indexlp) = vector2(indexlp)+vector1(mm)*vector1(nn)*valuelp
        !end if
        !if (indexlp1 /= 0) then
        valuelp1 = value_lpext1(iij)
        vector2(indexlp1) = vector2(indexlp1)+vector1(mm)*vector1(nn)*valuelp1
        !end if
        iij = iij+1
        mm = mm+1
      end do
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
      do i=1,lpext_wei(iii+3)
        indexlp = index_lpext(iij)
        indexlp1 = index_lpext1(iij)
        !if (indexlp /= 0) then
        valuelp = value_lpext(iij)
        vector2(indexlp) = vector2(indexlp)+vector1(mm)*vector1(nn)*valuelp
        !end if
        !if (indexlp1 /= 0) then
        valuelp1 = value_lpext1(iij)
        vector2(indexlp1) = vector2(indexlp1)+vector1(mm)*vector1(nn)*valuelp1
        !end if
        iij = iij+1
        mm = mm+1
      end do
    end do
  end if

end do

end subroutine inn_ext_ss_drl_loop_unpack_g

subroutine inn_ext_st_drl_loop_unpack_g(iplplwei,iplprwei)

use gugaci_global, only: ibsm_ext, idownwei_g131415, iesm_ext, index_lpext, index_lpext1, ism_g2g4, iwt_sm_s_ext, logic_g1415, &
                         logic_g2g4a, logic_g2g4b, logic_g34a, logic_g34b, logic_g35a, logic_g35b, logic_g36a, logic_g36b, &
                         lpend34a, lpend34b, lpend35a, lpend35b, lpend36a, lpend36b, lpext_wei, lpsta34a, lpsta34b, lpsta35a, &
                         lpsta35b, lpsta36a, lpsta36b, ng_sm, value_lpext, value_lpext1, vector1, vector2
use Symmetry_Info, only: Mul
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: iplplwei, iplprwei
integer(kind=iwp) :: i, ia, iaend, iasta, ib, ibend, ibsta, icle, ii, ii0, iii, iij, ilwtmp, indexlp, indexlp1, irwtmp, isma, &
                     ismb, lpend34, lpend35, lpend36, lpsta34, lpsta35, lpsta36, mm, mm0, nn, nn0, nna, nnb
real(kind=wp) :: valuelp, valuelp1

!write(u6,*) ' st_test 2/2'

ii = 1
if (logic_g1415) then
  mm = iplplwei
  nn = iplprwei
  do i=1,idownwei_g131415
    mm = mm+1
    nn = nn+1
    indexlp = index_lpext(ii)
    indexlp1 = index_lpext1(ii)
    !if (indexlp /= 0) then
    valuelp = value_lpext(ii)
    vector2(indexlp) = vector2(indexlp)+vector1(mm)*vector1(nn)*valuelp
    !end if
    !if (indexlp1 /= 0) then
    valuelp1 = value_lpext1(ii)
    vector2(indexlp1) = vector2(indexlp1)+vector1(mm)*vector1(nn)*valuelp1
    !end if
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

      do ia=iasta,min(iaend,ib-1)
        mm = mm+1
        nnb = nnb+1

        indexlp = index_lpext(ii)
        !if (indexlp /= 0) then
        valuelp = value_lpext(ii)
        vector2(indexlp) = vector2(indexlp)+vector1(mm)*vector1(nna)*valuelp
        !end if
        ii = ii+1

        indexlp = index_lpext(ii)
        !if (indexlp /= 0) then
        valuelp = value_lpext(ii)
        vector2(indexlp) = vector2(indexlp)+vector1(mm)*vector1(nnb)*valuelp
        !end if
        ii = ii+1
      end do
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
      do ia=iasta,min(iaend,ib-1)
        mm = mm+1
        nnb = nnb+1

        indexlp = index_lpext(ii)
        !if (indexlp /= 0) then
        valuelp = value_lpext(ii)
        vector2(indexlp) = vector2(indexlp)+vector1(mm)*vector1(nnb)*valuelp
        !end if
        ii = ii+1

        indexlp = index_lpext(ii)
        !if (indexlp /= 0) then
        valuelp = value_lpext(ii)
        vector2(indexlp) = vector2(indexlp)+vector1(mm)*vector1(nna)*valuelp
        !end if
        ii = ii+1

      end do
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
        indexlp = index_lpext(ii)
        !if (indexlp /= 0) then
        valuelp = value_lpext(ii)
        vector2(indexlp) = vector2(indexlp)+vector1(mm)*vector1(nn)*valuelp
        !end if
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

      do i=1,lpext_wei(iii+3)
        indexlp = index_lpext(iij)
        !if (indexlp /= 0) then
        valuelp = value_lpext(iij)
        vector2(indexlp) = vector2(indexlp)+vector1(mm)*vector1(nn)*valuelp
        !end if
        iij = iij+1
        mm = mm+1
      end do
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

      do i=1,lpext_wei(iii+3)
        indexlp = index_lpext(iij)
        !if (indexlp /= 0) then
        valuelp = -value_lpext(iij)
        vector2(indexlp) = vector2(indexlp)+vector1(mm)*vector1(nn)*valuelp
        !end if
        iij = iij+1
        mm = mm+1
      end do
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

      do i=1,lpext_wei(iii+3)
        indexlp = index_lpext(iij)
        !if (indexlp /= 0) then
        valuelp = -value_lpext(iij)
        vector2(indexlp) = vector2(indexlp)+vector1(mm)*vector1(nn)*valuelp
        !end if
        iij = iij+1
        mm = mm+1
      end do
    end do
  end if

end do

end subroutine inn_ext_st_drl_loop_unpack_g

subroutine inn_ext_tt_drl_loop_unpack_g(iplplwei,iplprwei,n1415)

use gugaci_global, only: idownwei_g131415, index_lpext, index_lpext1, logic_g1415, logic_g34a, logic_g34b, logic_g35a, logic_g35b, &
                         logic_g36a, logic_g36b, lpend34a, lpend34b, lpend35a, lpend35b, lpend36a, lpend36b, lpext_wei, lpsta34a, &
                         lpsta34b, lpsta35a, lpsta35b, lpsta36a, lpsta36b, value_lpext, value_lpext1, vector1, vector2
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: iplplwei, iplprwei, n1415
integer(kind=iwp) :: i, icle, ii, ii0, iii, iij, ilwtmp, indexlp, indexlp1, irwtmp, lpend34, lpend35, lpend36, lpsta34, lpsta35, &
                     lpsta36, mm, mm0, nn, nn0
real(kind=wp) :: valuelp, valuelp1
logical(kind=iwp) :: logic_g14150, logic_g34b0, logic_g35b0, logic_g36b0

!write(u6,*) '  tt_test 2/2'
logic_g14150 = logic_g1415
logic_g36b0 = logic_g36b
logic_g35b0 = logic_g35b
logic_g34b0 = logic_g34b
if (iplplwei == iplprwei) then
  logic_g14150 = .false.
  logic_g36b0 = .false.
  logic_g35b0 = .false.
  logic_g34b0 = .false.
end if
ii = 1
if (logic_g14150) then
  mm = iplplwei
  nn = iplprwei
  !ii = iista
  do i=1,idownwei_g131415
    mm = mm+1
    nn = nn+1
    indexlp = index_lpext(ii)
    !if (indexlp /= 0) then
    valuelp = value_lpext(ii)
    vector2(indexlp) = vector2(indexlp)+vector1(mm)*vector1(nn)*valuelp
    !end if
    indexlp1 = index_lpext1(ii)
    !if (indexlp1 /= 0) then
    valuelp1 = value_lpext1(ii)
    vector2(indexlp1) = vector2(indexlp1)+vector1(mm)*vector1(nn)*valuelp1
    !end if
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
        indexlp = index_lpext(ii)
        !if (indexlp /= 0) then
        valuelp = value_lpext(ii)
        vector2(indexlp) = vector2(indexlp)+vector1(mm)*vector1(nn)*valuelp
        !end if
        indexlp1 = index_lpext1(ii)
        !if (indexlp1 /= 0) then
        valuelp1 = value_lpext1(ii)
        vector2(indexlp1) = vector2(indexlp1)+vector1(mm)*vector1(nn)*valuelp1
        !end if
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
    else if ((icle == 2) .and. logic_g35b0) then
      mm0 = iplprwei
      nn0 = iplplwei
      lpsta35 = lpsta35b
      lpend35 = lpend35b
    end if
    do iii=lpsta35,lpend35,4
      mm = mm0+lpext_wei(iii)
      nn = nn0+lpext_wei(iii+1)
      iij = ii0+lpext_wei(iii+2)

      do i=1,lpext_wei(iii+3)
        indexlp = index_lpext(iij)
        !if (indexlp /= 0) then
        valuelp = -value_lpext(iij)
        vector2(indexlp) = vector2(indexlp)+vector1(mm)*vector1(nn)*valuelp
        !end if
        indexlp1 = index_lpext1(iij)
        !if (indexlp1 /= 0) then
        valuelp1 = -value_lpext1(iij)
        vector2(indexlp1) = vector2(indexlp1)+vector1(mm)*vector1(nn)*valuelp1
        !end if
        iij = iij+1
        mm = mm+1
      end do
    end do
  end if
  !cycle

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

      do i=1,lpext_wei(iii+3)
        indexlp = index_lpext(iij)
        !if (indexlp /= 0) then
        valuelp = value_lpext(iij)
        vector2(indexlp) = vector2(indexlp)+vector1(mm)*vector1(nn)*valuelp
        !end if
        indexlp1 = index_lpext1(iij)
        !if (indexlp1 /= 0) then
        valuelp1 = value_lpext1(iij)
        vector2(indexlp1) = vector2(indexlp1)+vector1(mm)*vector1(nn)*valuelp1
        !end if
        iij = iij+1
        mm = mm+1
      end do

    end do
  end if

end do

end subroutine inn_ext_tt_drl_loop_unpack_g

subroutine inn_ext_ts_drl_loop_unpack_g(iplplwei,iplprwei)

use gugaci_global, only: ibsm_ext, idownwei_g131415, iesm_ext, index_lpext, index_lpext1, ism_g2g4, iwt_sm_s_ext, logic_g1415, &
                         logic_g2g4a, logic_g2g4b, logic_g34a, logic_g34b, logic_g35a, logic_g35b, logic_g36a, logic_g36b, &
                         lpend34a, lpend34b, lpend35a, lpend35b, lpend36a, lpend36b, lpext_wei, lpsta34a, lpsta34b, lpsta35a, &
                         lpsta35b, lpsta36a, lpsta36b, ng_sm, value_lpext, value_lpext1, vector1, vector2
use Symmetry_Info, only: Mul
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: iplplwei, iplprwei
integer(kind=iwp) :: i, ia, iaend, iasta, ib, ibend, ibsta, icle, ii, ii0, iii, iij, ilwtmp, indexlp, indexlp1, irwtmp, isma, &
                     ismb, lpend34, lpend35, lpend36, lpsta34, lpsta35, lpsta36, mm, mm0, nn, nn0, nna, nnb
real(kind=wp) :: valuelp, valuelp1

!write(u6,*) '  ts_test 2/2'
ii = 1
!logic_g1415 = .false.
!logic_g2g4a = .false.
!logic_g2g4b = .false.
!logic_g36a = .false.
!logic_g36b = .false.
!logic_g35a = .false.
!logic_g35b = .false.
!logic_g34a = .false.
!logic_g34b = .false.
if (logic_g1415) then
  mm = iplplwei
  nn = iplprwei
  do i=1,idownwei_g131415
    mm = mm+1
    nn = nn+1
    indexlp = index_lpext(ii)
    !if (indexlp /= 0) then
    valuelp = value_lpext(ii)
    vector2(indexlp) = vector2(indexlp)+vector1(mm)*vector1(nn)*valuelp
    !end if
    indexlp1 = index_lpext1(ii)
    !if (indexlp1 /= 0) then
    valuelp1 = value_lpext1(ii)
    vector2(indexlp1) = vector2(indexlp1)+vector1(mm)*vector1(nn)*valuelp1
    !end if
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
      do ia=iasta,min(iaend,ib-1)
        mm = mm+1
        nnb = nnb+1
        indexlp = index_lpext(ii)
        !if (indexlp /= 0) then
        valuelp = value_lpext(ii)
        vector2(indexlp) = vector2(indexlp)+vector1(mm)*vector1(nna)*valuelp
        !end if
        ii = ii+1
        indexlp = index_lpext(ii)
        !if (indexlp /= 0) then
        valuelp = value_lpext(ii)
        vector2(indexlp) = vector2(indexlp)+vector1(mm)*vector1(nnb)*valuelp
        !end if
        ii = ii+1
      end do
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
      do ia=iasta,min(iaend,ib-1)
        mm = mm+1
        nnb = nnb+1
        indexlp = index_lpext(ii)
        !if (indexlp /= 0) then
        valuelp = value_lpext(ii)
        vector2(indexlp) = vector2(indexlp)+vector1(mm)*vector1(nnb)*valuelp
        !end if
        ii = ii+1
        indexlp = index_lpext(ii)
        !if (indexlp /= 0) then
        valuelp = value_lpext(ii)
        vector2(indexlp) = vector2(indexlp)+vector1(mm)*vector1(nna)*valuelp
        !end if
        ii = ii+1
      end do
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
        indexlp = index_lpext(ii)
        !if (indexlp /= 0) then
        valuelp = value_lpext(ii)
        vector2(indexlp) = vector2(indexlp)+vector1(mm)*vector1(nn)*valuelp
        !end if
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

      do i=1,lpext_wei(iii+3)
        indexlp = index_lpext(iij)
        !if (indexlp /= 0) then
        valuelp = -value_lpext(iij)
        vector2(indexlp) = vector2(indexlp)+vector1(mm)*vector1(nn)*valuelp
        !end if
        iij = iij+1
        mm = mm+1
      end do
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

      do i=1,lpext_wei(iii+3)
        indexlp = index_lpext(iij)
        !if (indexlp /= 0) then
        valuelp = value_lpext(iij)
        vector2(indexlp) = vector2(indexlp)+vector1(mm)*vector1(nn)*valuelp
        !end if
        iij = iij+1
        mm = mm+1
      end do
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

      do i=1,lpext_wei(iii+3)
        indexlp = index_lpext(iij)
        !if (indexlp /= 0) then
        valuelp = -value_lpext(iij)
        vector2(indexlp) = vector2(indexlp)+vector1(mm)*vector1(nn)*valuelp
        !end if
        iij = iij+1
        mm = mm+1
      end do
    end do
  end if
end do

end subroutine inn_ext_ts_drl_loop_unpack_g

subroutine inn_ext_ss_loop_unpack_g(iplplwei,iplprwei)

use gugaci_global, only: ibsm_ext, idownwei_g131415, iesm_ext, index_lpext, index_lpext1, ism_g2g4, iwt_sm_s_ext, logic_g1415, &
                         logic_g2g4a, logic_g2g4b, logic_g34a, logic_g34b, logic_g35a, logic_g35b, logic_g36a, logic_g36b, &
                         lpend34a, lpend34b, lpend35a, lpend35b, lpend36a, lpend36b, lpext_wei, lpsta34a, lpsta34b, lpsta35a, &
                         lpsta35b, lpsta36a, lpsta36b, ng_sm, nvalue_space_ss, value_lpext, value_lpext1, vector1, vector2
use Symmetry_Info, only: Mul
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: iplplwei, iplprwei
integer(kind=iwp) :: i, ia, iaend, iasta, ib, ibend, ibsta, icle, ii, ii0, iii, iij, ilwtmp, indexlp, indexlp1, irwtmp, isma, &
                     ismb, lpend34, lpend35, lpend36, lpsta34, lpsta35, lpsta36, mm, mm0, nn, nn0, nna, nnb
real(kind=wp) :: valuelp, valuelp1

!write(u6,*) ' ss_test 1/2'
ii = 1
if (logic_g1415) then
  mm = iplplwei
  nn = iplprwei
  do i=1,idownwei_g131415

    mm = mm+1
    nn = nn+1
    indexlp = index_lpext(ii)
    !if (indexlp /= 0) then
    valuelp = value_lpext(ii)
    vector2(indexlp) = vector2(indexlp)+vector1(mm)*vector1(nn)*valuelp
    !end if
    indexlp1 = index_lpext1(ii)
    !if (indexlp1 /= 0) then
    valuelp1 = value_lpext1(ii)
    vector2(indexlp1) = vector2(indexlp1)+vector1(mm)*vector1(nn)*valuelp1
    !end if
    ii = ii+1

    indexlp = index_lpext(ii)
    if (indexlp /= 0) then
      valuelp = value_lpext(ii)
      vector2(indexlp) = vector2(indexlp)+vector1(mm)*vector1(nn)*valuelp
    end if
    indexlp1 = index_lpext1(ii)
    if (indexlp1 /= 0) then
      valuelp1 = value_lpext1(ii)
      vector2(indexlp1) = vector2(indexlp1)+vector1(mm)*vector1(nn)*valuelp1
    end if
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
      do ia=iasta,min(iaend,ib-1)
        mm = mm+1
        nnb = nnb+1
        indexlp = index_lpext(ii)
        !if (indexlp /= 0) then
        valuelp = value_lpext(ii)
        vector2(indexlp) = vector2(indexlp)+vector1(mm)*vector1(nna)*valuelp
        !end if
        indexlp1 = index_lpext1(ii)
        !if (indexlp1 /= 0) then
        valuelp1 = value_lpext1(ii)
        vector2(indexlp1) = vector2(indexlp1)+vector1(mm)*vector1(nna)*valuelp1
        !end if
        ii = ii+1
        indexlp = index_lpext(ii)
        !if (indexlp /= 0) then
        valuelp = value_lpext(ii)
        vector2(indexlp) = vector2(indexlp)+vector1(mm)*vector1(nnb)*valuelp
        !end if
        indexlp1 = index_lpext1(ii)
        !if (indexlp1 /= 0) then
        valuelp1 = value_lpext1(ii)
        vector2(indexlp1) = vector2(indexlp1)+vector1(mm)*vector1(nnb)*valuelp1
        !end if
        ii = ii+1

      end do
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
      do ia=iasta,min(iaend,ib-1)
        mm = mm+1
        nnb = nnb+1
        indexlp = index_lpext(ii)
        !if (indexlp /= 0) then
        valuelp = value_lpext(ii)
        vector2(indexlp) = vector2(indexlp)+vector1(mm)*vector1(nnb)*valuelp
        !end if
        indexlp1 = index_lpext1(ii)
        !if (indexlp1 /= 0) then
        valuelp1 = value_lpext1(ii)
        vector2(indexlp1) = vector2(indexlp1)+vector1(mm)*vector1(nnb)*valuelp1
        !end if
        ii = ii+1
        indexlp = index_lpext(ii)
        !if (indexlp /= 0) then
        valuelp = value_lpext(ii)
        vector2(indexlp) = vector2(indexlp)+vector1(mm)*vector1(nna)*valuelp
        !end if
        indexlp1 = index_lpext1(ii)
        !if (indexlp1 /= 0) then
        valuelp1 = value_lpext1(ii)
        vector2(indexlp1) = vector2(indexlp1)+vector1(mm)*vector1(nna)*valuelp1
        !end if
        ii = ii+1

      end do
    end do
  end do
end if
!iaddii = ii-ii0)/2
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
        indexlp = index_lpext(ii)
        !if (indexlp /= 0) then
        valuelp = value_lpext(ii)
        vector2(indexlp) = vector2(indexlp)+vector1(mm)*vector1(nn)*valuelp
        !end if
        indexlp1 = index_lpext1(ii)
        !if (indexlp1 /= 0) then
        valuelp1 = value_lpext1(ii)
        vector2(indexlp1) = vector2(indexlp1)+vector1(mm)*vector1(nn)*valuelp1
        !end if

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
      do i=1,lpext_wei(iii+3)
        indexlp = index_lpext(iij)
        !if (indexlp /= 0) then
        valuelp = value_lpext(iij)
        vector2(indexlp) = vector2(indexlp)+vector1(mm)*vector1(nn)*valuelp
        !end if
        indexlp1 = index_lpext1(iij)
        !if (indexlp1 /= 0) then
        valuelp1 = value_lpext1(iij)
        vector2(indexlp1) = vector2(indexlp1)+vector1(mm)*vector1(nn)*valuelp1
        !end if
        iij = iij+1
        mm = mm+1
      end do
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
      do i=1,lpext_wei(iii+3)
        indexlp = index_lpext(iij)
        !if (indexlp /= 0) then
        valuelp = value_lpext(iij)
        vector2(indexlp) = vector2(indexlp)+vector1(mm)*vector1(nn)*valuelp
        !end if
        indexlp1 = index_lpext1(iij)
        !if (indexlp1 /= 0) then
        valuelp1 = value_lpext1(iij)
        vector2(indexlp1) = vector2(indexlp1)+vector1(mm)*vector1(nn)*valuelp1
        !end if
        iij = iij+1
        mm = mm+1
      end do
    end do
  end if
  ii0 = ii0+nvalue_space_ss
end do

end subroutine inn_ext_ss_loop_unpack_g

subroutine inn_ext_st_loop_unpack_g(iplplwei,iplprwei)

use gugaci_global, only: ibsm_ext, idownwei_g131415, iesm_ext, index_lpext, index_lpext1, ism_g2g4, iwt_sm_s_ext, logic_g1415, &
                         logic_g2g4a, logic_g2g4b, logic_g34a, logic_g34b, logic_g35a, logic_g35b, logic_g36a, logic_g36b, &
                         lpend34a, lpend34b, lpend35a, lpend35b, lpend36a, lpend36b, lpext_wei, lpsta34a, lpsta34b, lpsta35a, &
                         lpsta35b, lpsta36a, lpsta36b, ng_sm, nvalue_space_ss, value_lpext, value_lpext1, vector1, vector2
use Symmetry_Info, only: Mul
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: iplplwei, iplprwei
integer(kind=iwp) :: i, ia, iaend, iasta, ib, ibend, ibsta, icle, ii, ii0, iii, iij, ilwtmp, indexlp, indexlp1, irwtmp, isma, &
                     ismb, lpend34, lpend35, lpend36, lpsta34, lpsta35, lpsta36, mm, mm0, nn, nn0, nna, nnb
real(kind=wp) :: valuelp, valuelp1

!write(u6,*) ' st_test 1/2'
ii = 1
if (logic_g1415) then
  mm = iplplwei
  nn = iplprwei
  do i=1,idownwei_g131415
    mm = mm+1
    nn = nn+1
    indexlp = index_lpext(ii)
    !if (indexlp /= 0) then
    valuelp = value_lpext(ii)
    vector2(indexlp) = vector2(indexlp)+vector1(mm)*vector1(nn)*valuelp
    !end if
    indexlp1 = index_lpext1(ii)
    !if (indexlp1 /= 0) then
    valuelp1 = value_lpext1(ii)
    vector2(indexlp1) = vector2(indexlp1)+vector1(mm)*vector1(nn)*valuelp1
    !end if
    ii = ii+1

    indexlp = index_lpext(ii)
    if (indexlp /= 0) then
      valuelp = value_lpext(ii)
      vector2(indexlp) = vector2(indexlp)+vector1(mm)*vector1(nn)*valuelp
    end if
    indexlp1 = index_lpext1(ii)
    if (indexlp1 /= 0) then
      valuelp1 = value_lpext1(ii)
      vector2(indexlp1) = vector2(indexlp1)+vector1(mm)*vector1(nn)*valuelp1
    end if
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
      do ia=iasta,min(iaend,ib-1)
        mm = mm+1
        nnb = nnb+1
        indexlp = index_lpext(ii)
        !if (indexlp /= 0) then
        valuelp = value_lpext(ii)
        vector2(indexlp) = vector2(indexlp)+vector1(mm)*vector1(nna)*valuelp
        !end if
        indexlp1 = index_lpext1(ii)
        !if (indexlp1 /= 0) then
        valuelp1 = value_lpext1(ii)
        vector2(indexlp1) = vector2(indexlp1)+vector1(mm)*vector1(nna)*valuelp1
        !end if
        ii = ii+1

        indexlp = index_lpext(ii)
        !if (indexlp /= 0) then
        valuelp = value_lpext(ii)
        vector2(indexlp) = vector2(indexlp)+vector1(mm)*vector1(nnb)*valuelp
        !end if
        indexlp1 = index_lpext1(ii)
        !if (indexlp1 /= 0) then
        valuelp1 = value_lpext1(ii)
        vector2(indexlp1) = vector2(indexlp1)+vector1(mm)*vector1(nnb)*valuelp1
        !end if
        ii = ii+1
      end do
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
      do ia=iasta,min(iaend,ib-1)

        mm = mm+1
        nnb = nnb+1
        indexlp = index_lpext(ii)
        !if (indexlp /= 0) then
        valuelp = value_lpext(ii)
        vector2(indexlp) = vector2(indexlp)+vector1(mm)*vector1(nnb)*valuelp
        !end if
        indexlp1 = index_lpext1(ii)
        valuelp1 = value_lpext1(ii)
        !if (indexlp1 /= 0) then
        vector2(indexlp1) = vector2(indexlp1)+vector1(mm)*vector1(nnb)*valuelp1
        !end if
        ii = ii+1

        indexlp = index_lpext(ii)
        !if (indexlp /= 0) then
        valuelp = value_lpext(ii)
        vector2(indexlp) = vector2(indexlp)+vector1(mm)*vector1(nna)*valuelp
        !end if
        indexlp1 = index_lpext1(ii)
        !if (indexlp1 /= 0) then
        valuelp1 = value_lpext1(ii)
        vector2(indexlp1) = vector2(indexlp1)+vector1(mm)*vector1(nna)*valuelp1
        !end if
        ii = ii+1

      end do
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
        indexlp = index_lpext(ii)
        !if (indexlp /= 0) then
        valuelp = value_lpext(ii)
        vector2(indexlp) = vector2(indexlp)+vector1(mm)*vector1(nn)*valuelp
        !end if
        indexlp1 = index_lpext1(ii)
        !if (indexlp1 /= 0) then
        valuelp1 = value_lpext1(ii)
        vector2(indexlp1) = vector2(indexlp1)+vector1(mm)*vector1(nn)*valuelp1
        !end if
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
      do i=1,lpext_wei(iii+3)
        indexlp = index_lpext(iij)
        !if (indexlp /= 0) then
        valuelp = value_lpext(iij)
        vector2(indexlp) = vector2(indexlp)+vector1(mm)*vector1(nn)*valuelp
        !end if
        indexlp1 = index_lpext1(iij)
        !if (indexlp1 /= 0) then
        valuelp1 = value_lpext1(iij)
        vector2(indexlp1) = vector2(indexlp1)+vector1(mm)*vector1(nn)*valuelp1
        !end if
        iij = iij+1
        mm = mm+1
      end do
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
      do i=1,lpext_wei(iii+3)
        indexlp = index_lpext(iij)
        !if (indexlp /= 0) then
        valuelp = -value_lpext(iij)
        vector2(indexlp) = vector2(indexlp)+vector1(mm)*vector1(nn)*valuelp
        !end if
        indexlp1 = index_lpext1(iij)
        !if (indexlp1 /= 0) then
        valuelp1 = -value_lpext1(iij)
        vector2(indexlp1) = vector2(indexlp1)+vector1(mm)*vector1(nn)*valuelp1
        !end if
        iij = iij+1
        mm = mm+1
      end do
    end do
  end if

  !write(u6,*) 'st_g34',iplplwei,iplprwei,vector2(137)
  if (((icle == 1) .and. logic_g34a) .or. ((icle == 2) .and. logic_g34b)) then
    if (icle == 1) then
      mm0 = iplplwei
      nn0 = iplprwei
      lpsta34 = lpsta34a
      lpend34 = lpend34a
    else if ((icle == 2) .and. logic_g34b) then
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
      do i=1,lpext_wei(iii+3)
        indexlp = index_lpext(iij)
        !if (indexlp /= 0) then
        valuelp = -value_lpext(iij)
        vector2(indexlp) = vector2(indexlp)+vector1(mm)*vector1(nn)*valuelp
        !end if
        indexlp1 = index_lpext1(iij)
        !if (indexlp1 /= 0) then
        valuelp1 = -value_lpext1(iij)
        vector2(indexlp1) = vector2(indexlp1)+vector1(mm)*vector1(nn)*valuelp1
        !end if
        iij = iij+1
        mm = mm+1
      end do
    end do
  end if
  ii0 = ii0+nvalue_space_ss
end do

end subroutine inn_ext_st_loop_unpack_g

subroutine inn_ext_ts_loop_unpack_g(iplplwei,iplprwei)

use gugaci_global, only: ibsm_ext, idownwei_g131415, iesm_ext, index_lpext, index_lpext1, ism_g2g4, iwt_sm_s_ext, logic_g1415, &
                         logic_g2g4a, logic_g2g4b, logic_g34a, logic_g34b, logic_g35a, logic_g35b, logic_g36a, logic_g36b, &
                         lpend34a, lpend34b, lpend35a, lpend35b, lpend36a, lpend36b, lpext_wei, lpsta34a, lpsta34b, lpsta35a, &
                         lpsta35b, lpsta36a, lpsta36b, ng_sm, nvalue_space_ss, value_lpext, value_lpext1, vector1, vector2
use Symmetry_Info, only: Mul
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: iplplwei, iplprwei
integer(kind=iwp) :: i, ia, iaend, iasta, ib, ibend, ibsta, icle, ii, ii0, iii, iij, ilwtmp, indexlp, indexlp1, irwtmp, isma, &
                     ismb, lpend34, lpend35, lpend36, lpsta34, lpsta35, lpsta36, mm, mm0, nn, nn0, nna, nnb
real(kind=wp) :: valuelp, valuelp1

!write(u6,*) '  ts_test 1/2 '
ii = 1
if (logic_g1415) then
  mm = iplplwei
  nn = iplprwei
  do i=1,idownwei_g131415
    mm = mm+1
    nn = nn+1
    indexlp = index_lpext(ii)
    !if (indexlp /= 0) then
    valuelp = value_lpext(ii)
    vector2(indexlp) = vector2(indexlp)+vector1(mm)*vector1(nn)*valuelp
    !end if
    indexlp1 = index_lpext1(ii)
    !if (indexlp1 /= 0) then
    valuelp1 = value_lpext1(ii)
    vector2(indexlp1) = vector2(indexlp1)+vector1(mm)*vector1(nn)*valuelp1
    !end if
    ii = ii+1

    indexlp = index_lpext(ii)
    if (indexlp /= 0) then
      valuelp = value_lpext(ii)
      vector2(indexlp) = vector2(indexlp)+vector1(mm)*vector1(nn)*valuelp
    end if
    indexlp1 = index_lpext1(ii)
    if (indexlp1 /= 0) then
      valuelp1 = value_lpext1(ii)
      vector2(indexlp1) = vector2(indexlp1)+vector1(mm)*vector1(nn)*valuelp1
    end if
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
      do ia=iasta,min(iaend,ib-1)
        mm = mm+1
        nnb = nnb+1
        indexlp = index_lpext(ii)
        !if (indexlp /= 0) then
        valuelp = value_lpext(ii)
        vector2(indexlp) = vector2(indexlp)+vector1(mm)*vector1(nna)*valuelp
        !end if
        indexlp1 = index_lpext1(ii)
        !if (indexlp1 /= 0) then
        valuelp1 = value_lpext1(ii)
        vector2(indexlp1) = vector2(indexlp1)+vector1(mm)*vector1(nna)*valuelp1
        !end if
        ii = ii+1

        indexlp = index_lpext(ii)
        !if (indexlp /= 0) then
        valuelp = value_lpext(ii)
        vector2(indexlp) = vector2(indexlp)+vector1(mm)*vector1(nnb)*valuelp
        !end if
        indexlp1 = index_lpext1(ii)
        !if (indexlp1 /= 0) then
        valuelp1 = value_lpext1(ii)
        vector2(indexlp1) = vector2(indexlp1)+vector1(mm)*vector1(nnb)*valuelp1
        !end if
        ii = ii+1

      end do
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
      do ia=iasta,min(iaend,ib-1)
        mm = mm+1
        nnb = nnb+1
        indexlp = index_lpext(ii)
        !if (indexlp /= 0) then
        valuelp = value_lpext(ii)
        vector2(indexlp) = vector2(indexlp)+vector1(mm)*vector1(nnb)*valuelp
        !end if
        indexlp1 = index_lpext1(ii)
        !if (indexlp1 /= 0) then
        valuelp1 = value_lpext1(ii)
        vector2(indexlp1) = vector2(indexlp1)+vector1(mm)*vector1(nnb)*valuelp1
        !end if
        ii = ii+1

        indexlp = index_lpext(ii)
        !if (indexlp /= 0) then
        valuelp = value_lpext(ii)
        vector2(indexlp) = vector2(indexlp)+vector1(mm)*vector1(nna)*valuelp
        !end if
        indexlp1 = index_lpext1(ii)
        !if (indexlp1 /= 0) then
        valuelp1 = value_lpext1(ii)
        vector2(indexlp1) = vector2(indexlp1)+vector1(mm)*vector1(nna)*valuelp1
        !end if
        ii = ii+1

      end do
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
        indexlp = index_lpext(ii)
        !if (indexlp /= 0) then
        valuelp = value_lpext(ii)
        vector2(indexlp) = vector2(indexlp)+vector1(mm)*vector1(nn)*valuelp
        !end if
        indexlp1 = index_lpext1(ii)
        !if (indexlp1 /= 0) then
        valuelp1 = value_lpext1(ii)
        vector2(indexlp1) = vector2(indexlp1)+vector1(mm)*vector1(nn)*valuelp1
        !end if
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
      do i=1,lpext_wei(iii+3)
        indexlp = index_lpext(iij)
        !if (indexlp /= 0) then
        valuelp = -value_lpext(iij)
        vector2(indexlp) = vector2(indexlp)+vector1(mm)*vector1(nn)*valuelp
        !end if
        indexlp1 = index_lpext1(iij)
        !if (indexlp1 /= 0) then
        valuelp1 = -value_lpext1(iij)
        vector2(indexlp1) = vector2(indexlp1)+vector1(mm)*vector1(nn)*valuelp1
        !end if
        iij = iij+1
        mm = mm+1
      end do
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
      do i=1,lpext_wei(iii+3)
        indexlp = index_lpext(iij)
        !if (indexlp /= 0) then
        valuelp = value_lpext(iij)
        vector2(indexlp) = vector2(indexlp)+vector1(mm)*vector1(nn)*valuelp
        !end if
        indexlp1 = index_lpext1(iij)
        !if (indexlp1 /= 0) then
        valuelp1 = value_lpext1(iij)
        vector2(indexlp1) = vector2(indexlp1)+vector1(mm)*vector1(nn)*valuelp1
        !end if
        iij = iij+1
        mm = mm+1
      end do
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

      do i=1,lpext_wei(iii+3)
        indexlp = index_lpext(iij)
        !if (indexlp /= 0) then
        valuelp = -value_lpext(iij)
        vector2(indexlp) = vector2(indexlp)+vector1(mm)*vector1(nn)*valuelp
        !end if
        indexlp1 = index_lpext1(iij)
        !if (indexlp1 /= 0) then
        valuelp1 = -value_lpext1(iij)
        vector2(indexlp1) = vector2(indexlp1)+vector1(mm)*vector1(nn)*valuelp1
        !end if
        iij = iij+1
        mm = mm+1
      end do
    end do
  end if
  ii0 = ii0+nvalue_space_ss
end do

end subroutine inn_ext_ts_loop_unpack_g

subroutine inn_ext_tt_loop_unpack_g(iplplwei,iplprwei)

use gugaci_global, only: idownwei_g131415, index_lpext, index_lpext1, logic_g1415, logic_g34a, logic_g34b, logic_g35a, logic_g35b, &
                         logic_g36a, logic_g36b, lpend34a, lpend34b, lpend35a, lpend35b, lpend36a, lpend36b, lpext_wei, lpsta34a, &
                         lpsta34b, lpsta35a, lpsta35b, lpsta36a, lpsta36b, nvalue_space_ss, value_lpext, value_lpext1, vector1, &
                         vector2
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: iplplwei, iplprwei
integer(kind=iwp) :: i, icle, ii, ii0, iii, iij, ilwtmp, indexlp, indexlp1, irwtmp, lpend34, lpend35, lpend36, lpsta34, lpsta35, &
                     lpsta36, mm, mm0, nn, nn0
real(kind=wp) :: valuelp, valuelp1

!write(u6,*) '  tt_test 1/2'
ii = 1
if (logic_g1415) then
  mm = iplplwei
  nn = iplprwei

  do i=1,idownwei_g131415
    mm = mm+1
    nn = nn+1
    indexlp = index_lpext(ii)
    !if (indexlp /= 0) then
    valuelp = value_lpext(ii)
    vector2(indexlp) = vector2(indexlp)+vector1(mm)*vector1(nn)*valuelp
    !end if
    indexlp1 = index_lpext1(ii)
    !if (indexlp1 /= 0) then
    valuelp1 = value_lpext1(ii)
    vector2(indexlp1) = vector2(indexlp1)+vector1(mm)*vector1(nn)*valuelp1
    !end if

    ii = ii+1

    indexlp = index_lpext(ii)
    if (indexlp /= 0) then
      valuelp = value_lpext(ii)
      vector2(indexlp) = vector2(indexlp)+vector1(mm)*vector1(nn)*valuelp
    end if
    indexlp1 = index_lpext1(ii)
    if (indexlp1 /= 0) then
      valuelp1 = value_lpext1(ii)
      vector2(indexlp1) = vector2(indexlp1)+vector1(mm)*vector1(nn)*valuelp1
    end if
    ii = ii+1
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
        indexlp = index_lpext(ii)
        !if (indexlp /= 0) then
        valuelp = value_lpext(ii)
        vector2(indexlp) = vector2(indexlp)+vector1(mm)*vector1(nn)*valuelp
        !end if
        indexlp1 = index_lpext1(ii)
        !if (indexlp1 /= 0) then
        valuelp1 = value_lpext1(ii)
        vector2(indexlp1) = vector2(indexlp1)+vector1(mm)*vector1(nn)*valuelp1
        !end if
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

      do i=1,lpext_wei(iii+3)
        indexlp = index_lpext(iij)
        !if (indexlp /= 0) then
        valuelp = -value_lpext(iij)
        vector2(indexlp) = vector2(indexlp)+vector1(mm)*vector1(nn)*valuelp
        !end if
        indexlp1 = index_lpext1(iij)
        !if (indexlp1 /= 0) then
        valuelp1 = -value_lpext1(iij)
        vector2(indexlp1) = vector2(indexlp1)+vector1(mm)*vector1(nn)*valuelp1
        !end if

        iij = iij+1
        mm = mm+1
      end do
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

      do i=1,lpext_wei(iii+3)
        indexlp = index_lpext(iij)
        !if (indexlp /= 0) then
        valuelp = value_lpext(iij)
        vector2(indexlp) = vector2(indexlp)+vector1(mm)*vector1(nn)*valuelp
        !end if
        indexlp1 = index_lpext1(iij)
        !if (indexlp1 /= 0) then
        valuelp1 = value_lpext1(iij)
        vector2(indexlp1) = vector2(indexlp1)+vector1(mm)*vector1(nn)*valuelp1
        !end if

        iij = iij+1
        mm = mm+1
      end do
    end do
  end if
  ii0 = ii0+nvalue_space_ss
end do

end subroutine inn_ext_tt_loop_unpack_g

subroutine gsd_sequence_extspace_g(iplplwei,iplprwei)

use gugaci_global, only: index_lpext, index_lpext1, ivaluesta_g26, iweista_g25, iweista_g26, iweista_g28, logic_g25a, logic_g25b, &
                         logic_g26, logic_g28a, nint_g25, nint_g28, nwei_g25, nwei_g26, nwei_g28, v_sqtwo, value_lpext, &
                         value_lpext1, vector1, vector2
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: iplplwei, iplprwei
integer(kind=iwp) :: i, ilpvalue, indexlp, indexlp1, itmp, mm, nn, nn0
real(kind=wp) :: valuelp, valuelp1

!write(u6,*) '  sd_test 1/2','  ds_test 1'

ilpvalue = 0
if (logic_g25a) then
  mm = iplplwei+iweista_g25-1
  nn0 = iplprwei
  do itmp=1,nint_g25
    ilpvalue = ilpvalue+1
    indexlp = index_lpext(ilpvalue)
    valuelp = value_lpext(ilpvalue)
    indexlp1 = index_lpext1(ilpvalue)
    valuelp1 = value_lpext1(ilpvalue)
    nn = nn0
    do i=1,nwei_g25
      mm = mm+1
      nn = nn+1
      !if (indexlp /= 0) then
      vector2(indexlp) = vector2(indexlp)+vector1(mm)*vector1(nn)*valuelp
      !end if
      if (indexlp1 /= 0) then
        vector2(indexlp1) = vector2(indexlp1)+vector1(mm)*vector1(nn)*valuelp1
      end if
    end do
  end do
else if (logic_g25b) then
  mm = iplplwei+iweista_g25-1
  nn0 = iplprwei
  ilpvalue = ilpvalue+1
  do itmp=2,nint_g25
    ilpvalue = ilpvalue+1
    indexlp = index_lpext(ilpvalue)
    valuelp = value_lpext(ilpvalue)
    indexlp1 = index_lpext1(ilpvalue)
    valuelp1 = value_lpext1(ilpvalue)
    nn = nn0
    do i=1,itmp-1
      mm = mm+1
      nn = nn+1
      !if (indexlp /= 0) then
      vector2(indexlp) = vector2(indexlp)+vector1(mm)*vector1(nn)*valuelp
      !end if
      if (indexlp1 /= 0) then
        vector2(indexlp1) = vector2(indexlp1)+vector1(mm)*vector1(nn)*valuelp1
      end if
    end do
  end do

  mm = iplplwei+iweista_g28-1
  nn = iplprwei
  nn = nn+1
  do itmp=2,nwei_g28
    nn = nn+1
    ilpvalue = 0
    do i=1,itmp-1
      ilpvalue = ilpvalue+1
      mm = mm+1
      indexlp = index_lpext(ilpvalue)
      !if (indexlp /= 0) then
      valuelp = value_lpext(ilpvalue)
      vector2(indexlp) = vector2(indexlp)+vector1(mm)*vector1(nn)*valuelp
      !end if
      indexlp1 = index_lpext1(ilpvalue)
      if (indexlp1 /= 0) then
        valuelp1 = value_lpext1(ilpvalue)
        vector2(indexlp1) = vector2(indexlp1)+vector1(mm)*vector1(nn)*valuelp1
      end if
    end do
  end do
else if (logic_g28a) then
  mm = iplplwei+iweista_g28-1
  nn0 = iplprwei
  do nn=nn0+1,nn0+nwei_g28
    ilpvalue = 0
    do i=1,nint_g28
      ilpvalue = ilpvalue+1
      mm = mm+1
      indexlp = index_lpext(ilpvalue)
      !if (indexlp /= 0) then
      valuelp = value_lpext(ilpvalue)
      vector2(indexlp) = vector2(indexlp)+vector1(mm)*vector1(nn)*valuelp
      !end if
      indexlp1 = index_lpext1(ilpvalue)
      if (indexlp1 /= 0) then
        valuelp1 = value_lpext1(ilpvalue)
        vector2(indexlp1) = vector2(indexlp1)+vector1(mm)*vector1(nn)*valuelp1
      end if
    end do
  end do
end if

if (logic_g26) then
  ilpvalue = ivaluesta_g26
  mm = iplplwei+iweista_g26
  nn0 = iplprwei
  do nn=nn0+1,nn0+nwei_g26
    ilpvalue = ilpvalue+1
    indexlp = index_lpext(ilpvalue)
    !if (indexlp /= 0) then
    valuelp = value_lpext(ilpvalue)*v_sqtwo
    vector2(indexlp) = vector2(indexlp)+vector1(mm)*vector1(nn)*valuelp
    !end if
    indexlp1 = index_lpext1(ilpvalue)
    if (indexlp1 /= 0) then
      valuelp1 = value_lpext1(ilpvalue)*v_sqtwo
      vector2(indexlp1) = vector2(indexlp1)+vector1(mm)*vector1(nn)*valuelp1
    end if
    mm = mm+1
  end do
end if
!write(u6,*)  'out sd 1'

end subroutine gsd_sequence_extspace_g

subroutine gtd_sequence_extspace_g(iplplwei,iplprwei)

use gugaci_global, only: index_lpext, index_lpext1, iweista_g25, iweista_g28, logic_g25a, logic_g25b, logic_g28a, nint_g25, &
                         nint_g28, nwei_g25, nwei_g28, value_lpext, value_lpext1, vector1, vector2
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: iplplwei, iplprwei
integer(kind=iwp) :: i, ilpvalue, indexlp, indexlp1, itmp, mm, nn, nn0
real(kind=wp) :: valuelp, valuelp1

!write(u6,*) ' td_test _1/2',' dt_test '
ilpvalue = 0
if (logic_g25a) then
  mm = iplplwei+iweista_g25-1
  nn0 = iplprwei
  do itmp=1,nint_g25
    ilpvalue = ilpvalue+1
    indexlp = index_lpext(ilpvalue)
    valuelp = value_lpext(ilpvalue)
    indexlp1 = index_lpext1(ilpvalue)
    valuelp1 = value_lpext1(ilpvalue)

    nn = nn0
    do i=1,nwei_g25
      mm = mm+1
      nn = nn+1
      !if (indexlp /= 0) then
      vector2(indexlp) = vector2(indexlp)+vector1(mm)*vector1(nn)*valuelp
      !end if
      if (indexlp1 /= 0) then
        vector2(indexlp1) = vector2(indexlp1)+vector1(mm)*vector1(nn)*valuelp1
      end if
    end do
  end do
else if (logic_g25b) then
  mm = iplplwei+iweista_g25-1
  nn0 = iplprwei
  ilpvalue = ilpvalue+1
  do itmp=2,nint_g25
    ilpvalue = ilpvalue+1
    indexlp = index_lpext(ilpvalue)
    valuelp = value_lpext(ilpvalue)
    indexlp1 = index_lpext1(ilpvalue)
    valuelp1 = value_lpext1(ilpvalue)

    nn = nn0
    do i=1,itmp-1
      mm = mm+1
      nn = nn+1
      !if (indexlp /= 0) then
      vector2(indexlp) = vector2(indexlp)+vector1(mm)*vector1(nn)*valuelp
      !end if
      if (indexlp1 /= 0) then
        vector2(indexlp1) = vector2(indexlp1)+vector1(mm)*vector1(nn)*valuelp1
      end if
    end do
  end do
  mm = iplplwei+iweista_g28-1
  nn = iplprwei
  nn = nn+1
  do itmp=2,nwei_g28
    nn = nn+1
    ilpvalue = 0
    do i=1,itmp-1
      ilpvalue = ilpvalue+1
      mm = mm+1
      indexlp = index_lpext(ilpvalue)
      !if (indexlp /= 0) then
      valuelp = -value_lpext(ilpvalue)
      vector2(indexlp) = vector2(indexlp)+vector1(mm)*vector1(nn)*valuelp
      !end if
      indexlp1 = index_lpext1(ilpvalue)
      if (indexlp1 /= 0) then
        valuelp1 = -value_lpext1(ilpvalue)
        vector2(indexlp1) = vector2(indexlp1)+vector1(mm)*vector1(nn)*valuelp1
      end if
    end do
  end do
else if (logic_g28a) then
  mm = iplplwei+iweista_g28-1
  nn0 = iplprwei
  do nn=nn0+1,nn0+nwei_g28

    ilpvalue = 0
    do i=1,nint_g28
      ilpvalue = ilpvalue+1
      mm = mm+1
      indexlp = index_lpext(ilpvalue)
      !if (indexlp /= 0) then
      valuelp = -value_lpext(ilpvalue)
      vector2(indexlp) = vector2(indexlp)+vector1(mm)*vector1(nn)*valuelp
      !end if
      indexlp1 = index_lpext1(ilpvalue)
      if (indexlp1 /= 0) then
        valuelp1 = -value_lpext1(ilpvalue)
        vector2(indexlp1) = vector2(indexlp1)+vector1(mm)*vector1(nn)*valuelp1
      end if
    end do
  end do
end if

end subroutine gtd_sequence_extspace_g

subroutine gdv_sequence_extspace_g(ilw,irw)

use gugaci_global, only: ilsegdownwei, index_lpext, index_lpext1, value_lpext, value_lpext1, vector1, vector2
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: ilw, irw
integer(kind=iwp) :: iij, indexlp, indexlp1, mm, nn
real(kind=wp) :: valuelp, valuelp1

mm = ilw
nn = irw+1

do iij=1,ilsegdownwei
  mm = mm+1
  indexlp = index_lpext(iij)
  !if (indexlp /= 0) then
  valuelp = value_lpext(iij)
  vector2(indexlp) = vector2(indexlp)+vector1(mm)*vector1(nn)*valuelp
  !end if
  indexlp1 = index_lpext1(iij)
  if (indexlp1 /= 0) then
    valuelp1 = value_lpext1(iij)
    vector2(indexlp1) = vector2(indexlp1)+vector1(mm)*vector1(nn)*valuelp1
  end if
end do

end subroutine gdv_sequence_extspace_g

subroutine complete_sd_ar_ext_loop_g(ilwei,irwei,isdownwei)

use gugaci_global, only: icano_nnend, icano_nnsta, index_lpext, index_lpext1, value_lpext, value_lpext1, vector1, vector2
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: ilwei, irwei, isdownwei
integer(kind=iwp) :: ilpvalue, indexlp, indexlp1, mm, mm0, mmtmp, nn, nntmp
real(kind=wp) :: valuelp, valuelp1

!write(u6,*) 'sd_test 2/2','  td_test_2/2 2'
ilpvalue = 0
mm0 = ilwei
nn = irwei+icano_nnsta-1
do nntmp=icano_nnsta,icano_nnend
  nn = nn+1
  mm = mm0
  do mmtmp=1,isdownwei
    ilpvalue = ilpvalue+1
    mm = mm+1
    indexlp = index_lpext(ilpvalue)
    valuelp = value_lpext(ilpvalue)
    vector2(indexlp) = vector2(indexlp)+vector1(mm)*vector1(nn)*valuelp
    indexlp1 = index_lpext1(ilpvalue)
    if (indexlp1 /= 0) then

      valuelp1 = value_lpext1(ilpvalue)
      vector2(indexlp1) = vector2(indexlp1)+vector1(mm)*vector1(nn)*valuelp1
    end if
  end do
end do
!write(u6,*) 'out sd 2'

return

end subroutine complete_sd_ar_ext_loop_g

subroutine gdv_sequence_extspace1_g(ilw,irw,n)

use gugaci_global, only: dm1tmp, ilsegdownwei, index_lpext3, index_lpext4, index_lpext5, value_lpext3, value_lpext4, value_lpext5, &
                         vector1, vector2
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: ilw, irw, n
integer(kind=iwp) :: ilpvalue, indexlp2, indexlp3, indexlp4, mm, nii, nn
real(kind=wp) :: valuelp2, valuelp3, valuelp4

mm = ilw
nn = irw+1

do ilpvalue=1,ilsegdownwei
  mm = mm+1
  indexlp2 = index_lpext5(ilpvalue)
  valuelp2 = value_lpext5(ilpvalue)
  dm1tmp(indexlp2) = dm1tmp(indexlp2)+vector1(mm)*vector1(nn)*valuelp2

  do nii=1,n
    indexlp3 = index_lpext3(ilpvalue,nii)
    !if (indexlp3 /= 0) then
    valuelp3 = value_lpext3(ilpvalue,nii)
    vector2(indexlp3) = vector2(indexlp3)+vector1(mm)*vector1(nn)*valuelp3
    !end if

    indexlp4 = index_lpext4(ilpvalue,nii)
    if (indexlp4 /= 0) then
      valuelp4 = value_lpext4(ilpvalue,nii)
      vector2(indexlp4) = vector2(indexlp4)+vector1(mm)*vector1(nn)*valuelp4
    end if
  end do
end do

end subroutine gdv_sequence_extspace1_g

subroutine gtd_sequence_extspace1_g(iplplwei,iplprwei,n)

use gugaci_global, only: dm1tmp, index_lpext3, index_lpext4, index_lpext5, iweista_g25, iweista_g28, logic_g25a, logic_g25b, &
                         logic_g28a, nint_g25, nint_g28, nwei_g25, nwei_g28, value_lpext3, value_lpext4, value_lpext5, vector1, &
                         vector2
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: iplplwei, iplprwei, n
integer(kind=iwp) :: i, ilpvalue, indexlp2, indexlp3, indexlp4, itmp, mm, nii, nn, nn0
real(kind=wp) :: valuelp2, valuelp3, valuelp4

!write(u6,*) ' td_test _1/2',' dt_test '
ilpvalue = 0
if (logic_g25a) then
  mm = iplplwei+iweista_g25-1
  nn0 = iplprwei
  do itmp=1,nint_g25
    ilpvalue = ilpvalue+1
    indexlp2 = index_lpext5(ilpvalue)
    valuelp2 = value_lpext5(ilpvalue)
    nn = nn0
    do i=1,nwei_g25
      mm = mm+1
      nn = nn+1
      dm1tmp(indexlp2) = dm1tmp(indexlp2)+vector1(mm)*vector1(nn)*valuelp2
      do nii=1,n
        indexlp3 = index_lpext3(ilpvalue,nii)
        !if (indexlp3 /= 0) then
        valuelp3 = value_lpext3(ilpvalue,nii)
        vector2(indexlp3) = vector2(indexlp3)+vector1(mm)*vector1(nn)*valuelp3
        !end if
        indexlp4 = index_lpext4(ilpvalue,nii)
        if (indexlp4 /= 0) then
          valuelp4 = value_lpext4(ilpvalue,nii)
          vector2(indexlp4) = vector2(indexlp4)+vector1(mm)*vector1(nn)*valuelp4
        end if
      end do
    end do
  end do
else if (logic_g25b) then
  mm = iplplwei+iweista_g25-1
  nn0 = iplprwei
  ilpvalue = ilpvalue+1
  do itmp=2,nint_g25
    ilpvalue = ilpvalue+1
    indexlp2 = index_lpext5(ilpvalue)
    valuelp2 = value_lpext5(ilpvalue)

    nn = nn0
    do i=1,itmp-1
      mm = mm+1
      nn = nn+1
      dm1tmp(indexlp2) = dm1tmp(indexlp2)+vector1(mm)*vector1(nn)*valuelp2
      do nii=1,n
        indexlp3 = index_lpext3(ilpvalue,nii)
        !if (indexlp3 /= 0) then
        valuelp3 = value_lpext3(ilpvalue,nii)
        vector2(indexlp3) = vector2(indexlp3)+vector1(mm)*vector1(nn)*valuelp3
        !end if
        indexlp4 = index_lpext4(ilpvalue,nii)
        if (indexlp4 /= 0) then
          valuelp4 = value_lpext4(ilpvalue,nii)
          vector2(indexlp4) = vector2(indexlp4)+vector1(mm)*vector1(nn)*valuelp4
        end if
      end do
    end do
  end do
  mm = iplplwei+iweista_g28-1
  nn = iplprwei
  nn = nn+1
  do itmp=2,nwei_g28
    nn = nn+1
    ilpvalue = 0
    do i=1,itmp-1
      ilpvalue = ilpvalue+1
      indexlp2 = index_lpext5(ilpvalue)
      valuelp2 = -value_lpext5(ilpvalue)
      mm = mm+1
      dm1tmp(indexlp2) = dm1tmp(indexlp2)+vector1(mm)*vector1(nn)*valuelp2
      do nii=1,n
        indexlp3 = index_lpext3(ilpvalue,nii)
        !if (indexlp3 /= 0) then
        valuelp3 = -value_lpext3(ilpvalue,nii)
        vector2(indexlp3) = vector2(indexlp3)+vector1(mm)*vector1(nn)*valuelp3
        !end if
        indexlp4 = index_lpext4(ilpvalue,nii)
        if (indexlp4 /= 0) then
          valuelp4 = -value_lpext4(ilpvalue,nii)
          vector2(indexlp4) = vector2(indexlp4)+vector1(mm)*vector1(nn)*valuelp4
        end if
      end do
    end do
  end do
else if (logic_g28a) then
  mm = iplplwei+iweista_g28-1
  nn0 = iplprwei
  do nn=nn0+1,nn0+nwei_g28

    ilpvalue = 0
    do i=1,nint_g28
      ilpvalue = ilpvalue+1
      indexlp2 = index_lpext5(ilpvalue)
      valuelp2 = -value_lpext5(ilpvalue)
      mm = mm+1
      dm1tmp(indexlp2) = dm1tmp(indexlp2)+vector1(mm)*vector1(nn)*valuelp2
      do nii=1,n
        indexlp3 = index_lpext3(ilpvalue,nii)
        !if (indexlp3 /= 0) then
        valuelp3 = -value_lpext3(ilpvalue,nii)
        vector2(indexlp3) = vector2(indexlp3)+vector1(mm)*vector1(nn)*valuelp3
        !end if
        indexlp4 = index_lpext4(ilpvalue,nii)
        if (indexlp4 /= 0) then
          valuelp4 = -value_lpext4(ilpvalue,nii)
          vector2(indexlp4) = vector2(indexlp4)+vector1(mm)*vector1(nn)*valuelp4
        end if
      end do
    end do
  end do
end if

end subroutine gtd_sequence_extspace1_g

subroutine gsd_sequence_extspace1_g(iplplwei,iplprwei,n)

use gugaci_global, only: dm1tmp, index_lpext3, index_lpext4, index_lpext5, ivaluesta_g26, iweista_g25, iweista_g26, iweista_g28, &
                         logic_g25a, logic_g25b, logic_g26, logic_g28a, nint_g25, nint_g28, nwei_g25, nwei_g26, nwei_g28, &
                         v_sqtwo, value_lpext3, value_lpext4, value_lpext5, vector1, vector2
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: iplplwei, iplprwei, n
integer(kind=iwp) :: i, ilpvalue, indexlp2, indexlp3, indexlp4, itmp, mm, nii, nn, nn0
real(kind=wp) :: valuelp2, valuelp3, valuelp4

!write(u6,*) '  sd_test 1/2','  ds_test 0'

ilpvalue = 0
if (logic_g25a) then
  mm = iplplwei+iweista_g25-1
  nn0 = iplprwei
  do itmp=1,nint_g25
    ilpvalue = ilpvalue+1
    indexlp2 = index_lpext5(ilpvalue)
    valuelp2 = value_lpext5(ilpvalue)

    nn = nn0
    do i=1,nwei_g25
      mm = mm+1
      nn = nn+1
      dm1tmp(indexlp2) = dm1tmp(indexlp2)+vector1(mm)*vector1(nn)*valuelp2
      do nii=1,n
        indexlp3 = index_lpext3(ilpvalue,nii)
        !if (indexlp3 /= 0) then
        valuelp3 = value_lpext3(ilpvalue,nii)
        vector2(indexlp3) = vector2(indexlp3)+vector1(mm)*vector1(nn)*valuelp3
        !end if
        indexlp4 = index_lpext4(ilpvalue,nii)
        if (indexlp4 /= 0) then
          valuelp4 = value_lpext4(ilpvalue,nii)
          vector2(indexlp4) = vector2(indexlp4)+vector1(mm)*vector1(nn)*valuelp4
        end if
      end do
    end do
  end do
else if (logic_g25b) then
  mm = iplplwei+iweista_g25-1
  nn0 = iplprwei
  ilpvalue = ilpvalue+1
  do itmp=2,nint_g25
    ilpvalue = ilpvalue+1
    indexlp2 = index_lpext5(ilpvalue)
    valuelp2 = value_lpext5(ilpvalue)

    nn = nn0
    do i=1,itmp-1
      mm = mm+1
      nn = nn+1
      dm1tmp(indexlp2) = dm1tmp(indexlp2)+vector1(mm)*vector1(nn)*valuelp2
      do nii=1,n
        indexlp3 = index_lpext3(ilpvalue,nii)
        !if (indexlp3 /= 0) then
        valuelp3 = value_lpext3(ilpvalue,nii)
        vector2(indexlp3) = vector2(indexlp3)+vector1(mm)*vector1(nn)*valuelp3
        !end if
        indexlp4 = index_lpext4(ilpvalue,nii)
        if (indexlp4 /= 0) then
          valuelp4 = value_lpext4(ilpvalue,nii)
          vector2(indexlp4) = vector2(indexlp4)+vector1(mm)*vector1(nn)*valuelp4
        end if
      end do

    end do
  end do

  mm = iplplwei+iweista_g28-1
  nn = iplprwei
  nn = nn+1
  do itmp=2,nwei_g28
    nn = nn+1
    ilpvalue = 0
    do i=1,itmp-1
      ilpvalue = ilpvalue+1
      indexlp2 = index_lpext5(ilpvalue)
      valuelp2 = value_lpext5(ilpvalue)
      mm = mm+1
      dm1tmp(indexlp2) = dm1tmp(indexlp2)+vector1(mm)*vector1(nn)*valuelp2
      do nii=1,n
        indexlp3 = index_lpext3(ilpvalue,nii)
        !if (indexlp3 /= 0) then
        valuelp3 = value_lpext3(ilpvalue,nii)
        vector2(indexlp3) = vector2(indexlp3)+vector1(mm)*vector1(nn)*valuelp3
        !end if
        indexlp4 = index_lpext4(ilpvalue,nii)
        if (indexlp4 /= 0) then
          valuelp4 = value_lpext4(ilpvalue,nii)
          vector2(indexlp4) = vector2(indexlp4)+vector1(mm)*vector1(nn)*valuelp4
        end if
      end do
    end do
  end do
else if (logic_g28a) then
  mm = iplplwei+iweista_g28-1
  nn0 = iplprwei
  do nn=nn0+1,nn0+nwei_g28
    ilpvalue = 0
    do i=1,nint_g28
      ilpvalue = ilpvalue+1
      indexlp2 = index_lpext5(ilpvalue)
      valuelp2 = value_lpext5(ilpvalue)
      mm = mm+1
      dm1tmp(indexlp2) = dm1tmp(indexlp2)+vector1(mm)*vector1(nn)*valuelp2
      do nii=1,n
        indexlp3 = index_lpext3(ilpvalue,nii)
        !if (indexlp3 /= 0) then
        valuelp3 = value_lpext3(ilpvalue,nii)
        vector2(indexlp3) = vector2(indexlp3)+vector1(mm)*vector1(nn)*valuelp3
        !end if
        indexlp4 = index_lpext4(ilpvalue,nii)
        if (indexlp4 /= 0) then
          valuelp4 = value_lpext4(ilpvalue,nii)
          vector2(indexlp4) = vector2(indexlp4)+vector1(mm)*vector1(nn)*valuelp4
        end if
      end do
    end do
  end do
end if

if (logic_g26) then
  ilpvalue = ivaluesta_g26
  mm = iplplwei+iweista_g26
  nn0 = iplprwei
  do nn=nn0+1,nn0+nwei_g26
    ilpvalue = ilpvalue+1
    indexlp2 = index_lpext5(ilpvalue)
    valuelp2 = value_lpext5(ilpvalue)*v_sqtwo

    dm1tmp(indexlp2) = dm1tmp(indexlp2)+vector1(mm)*vector1(nn)*valuelp2
    do nii=1,n
      indexlp3 = index_lpext3(ilpvalue,nii)
      !if (indexlp3 /= 0) then
      valuelp3 = value_lpext3(ilpvalue,nii)*v_sqtwo
      vector2(indexlp3) = vector2(indexlp3)+vector1(mm)*vector1(nn)*valuelp3
      !end if
      indexlp4 = index_lpext4(ilpvalue,nii)
      if (indexlp4 /= 0) then
        valuelp4 = value_lpext4(ilpvalue,nii)*v_sqtwo
        vector2(indexlp4) = vector2(indexlp4)+vector1(mm)*vector1(nn)*valuelp4
      end if
    end do
    mm = mm+1
  end do
end if
!write(u6,*) 'out ds 0'

end subroutine gsd_sequence_extspace1_g

subroutine inn_ext_dd_loop_unpack_g(iplplwei,iplprwei)

use gugaci_global, only: ildownwei_segdd, index_lpext, index_lpext1, int_dd_drl, irdownwei_segdd, logic_g49a, logic_g49b, &
                         logic_g50, value_lpext, value_lpext1, vector1, vector2
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: iplplwei, iplprwei
integer(kind=iwp) :: i, icle, ii, ildownwei, indexlp, indexlp1, irdownwei, j, mm, mm0, nn
real(kind=wp) :: valuelp, valuelp1

!write(nf2,*) 'logic_g49b',logic_g50,logic_g49a,logic_g49b

ii = 1
if (logic_g50) then
  if (logic_g49b) then
    mm = iplplwei
    nn = iplprwei
    do i=1,ildownwei_segdd
      mm = mm+1
      nn = nn+1
      indexlp = index_lpext(ii)
      if (indexlp /= 0) then
        valuelp = value_lpext(ii)
        vector2(indexlp) = vector2(indexlp)+vector1(mm)*vector1(nn)*valuelp
      end if
      indexlp1 = index_lpext1(ii)
      if (indexlp1 /= 0) then
        valuelp1 = value_lpext1(ii)
        vector2(indexlp1) = vector2(indexlp1)+vector1(mm)*vector1(nn)*valuelp1
      end if
      ii = ii+1
    end do
  end if

  ii = ii+int_dd_drl
  mm0 = iplplwei
  nn = iplprwei+1
  do icle=1,2
    do j=2,ildownwei_segdd
      nn = nn+1
      mm = mm0
      do i=1,j-1
        mm = mm+1
        indexlp = index_lpext(ii)
        if (indexlp /= 0) then
          valuelp = value_lpext(ii)
          vector2(indexlp) = vector2(indexlp)+vector1(mm)*vector1(nn)*valuelp
        end if
        indexlp1 = index_lpext1(ii)
        if (indexlp1 /= 0) then
          valuelp1 = value_lpext1(ii)
          vector2(indexlp1) = vector2(indexlp1)+vector1(mm)*vector1(nn)*valuelp1
        end if
        ii = ii+1
      end do
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
    mm = mm0
    do i=1,ildownwei
      mm = mm+1
      indexlp = index_lpext(ii)
      !if (indexlp /= 0) then
      valuelp = value_lpext(ii)
      vector2(indexlp) = vector2(indexlp)+vector1(mm)*vector1(nn)*valuelp
      !end if
      indexlp1 = index_lpext1(ii)
      !if (indexlp1 /= 0) then
      valuelp1 = value_lpext1(ii)
      vector2(indexlp1) = vector2(indexlp1)+vector1(mm)*vector1(nn)*valuelp1
      !end if
      ii = ii+1
    end do
  end do
end if

end subroutine inn_ext_dd_loop_unpack_g

subroutine inn_ext_sv_loop_unpack_g(ilw,irw)

use gugaci_global, only: ilsegdownwei, index_lpext, index_lpext1, value_lpext, value_lpext1, vector1, vector2
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: ilw, irw
integer(kind=iwp) :: iij, indexlp, indexlp1, mm, nn
real(kind=wp) :: valuelp, valuelp1

mm = ilw
nn = irw+1
do iij=1,ilsegdownwei
  mm = mm+1
  indexlp = index_lpext(iij)
  !if (indexlp /= 0) then
  valuelp = value_lpext(iij)
  vector2(indexlp) = vector2(indexlp)+vector1(mm)*vector1(nn)*valuelp
  !end if
  indexlp1 = index_lpext1(iij)
  if (indexlp1 /= 0) then
    valuelp1 = value_lpext1(iij)
    vector2(indexlp1) = vector2(indexlp1)+vector1(mm)*vector1(nn)*valuelp1
  end if
end do

end subroutine inn_ext_sv_loop_unpack_g

! density matrix formation : vector2=ci*<i|epq,rs|j>*cj and dm1=ci*<i|e
! for norb_act<>0         mg1,mg2,mg3,mg4,mg5:
! idb=1  in dbl_space      ity_up=0-5               0 ,jpad,iwdl,iwdr,
! idb=2  in act_space      ity_up=0-5,itdown=0,3      jph, jpe,iwal,iwa
! idb=3  between dbl and act   ity_up=0-5,itdown=0,3      jpe,iwdl,iwdr,

! this subroutine prodab_1 does the dm1 part, which corresponds to voin
subroutine prodab_1(idb,mg1,mg2,mg3,mg4,mg5,wl,jpr,mg6,mg7)

use gugaci_global, only: log_prod
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: idb, mg1, mg2, mg3, mg4, mg5, jpr, mg6, mg7
real(kind=wp), intent(in) :: wl

if (log_prod == 1) call prodab_h_1(idb,mg1,mg2,mg3,mg4,mg5,wl,jpr,mg6,mg7)
if (log_prod == 2) call prodab_h0_1(idb,mg1,mg2,mg3,mg4,mg5,wl,jpr,mg6,mg7)

return

end subroutine prodab_1

subroutine prodab_h_1(idb,mg1,mg2,mg3,mg4,mg5,wl,jpr,mg6,mg7)

use gugaci_global, only: dm1tmp, ican_a, ihy, ipae, ipael, iseg_downwei, iw_downwei, iy, jpad, jpad_upwei, jpadl, jphy, loputmp, &
                         nu_ae, vector1
use stdalloc, only: mma_allocate, mma_deallocate
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: idb, mg1, mg2, mg3, mg4, mg5, jpr, mg6, mg7
real(kind=wp), intent(in) :: wl
integer(kind=iwp) :: ii, in_, ipae_, isegdownwei, iwa, iwadl, iwadr, iwal, iwar, iwd, iwdl, iwdown, iwdr, iwe, iwl, iwr, iwupwei, &
                     jpe, jph, jpl, jpy, jwd, jwnu, jwu, lp, lwnu, m, mg67, mm, mp, nn
integer(kind=iwp), allocatable :: lopu(:,:)
integer(kind=iwp), external :: iwalk_ad

select case (idb)
  case default ! (1)
    ! in dbl_space
    jpad = mg2
    iwdl = mg3
    iwdr = mg4
    do ipae_=1,25
      ipae = ipae_ ! ipae is in global module, is this necessary?
      if (nu_ae(ipae) == 0) cycle
      iwdown = iw_downwei(jpad,ipae)
      if (iwdown == 0) cycle
      lwnu = iseg_downwei(ipae)
      do iwa=0,iwdown-1
        iwadl = iwalk_ad(jpad,ipae,iwa,iwdl)
        iwadr = iwalk_ad(jpad,ipae,iwa,iwdr)
        mm = iwadl
        nn = iwadr
        do m=1,lwnu
          mm = mm+1
          nn = nn+1
          !vector2(mm) = vector2(mm)+vector1(nn)*wl
          !vector2(nn) = vector2(nn)+vector1(mm)*wl
          !vector2(mm) = vector2(mm)+vector1(nn)*wl
          !vector2(nn) = vector2(nn)+vector1(mm)*wl
          mg67 = ican_a(mg7)+mg6
          dm1tmp(mg67) = dm1tmp(mg67)+vector1(nn)*wl*vector1(mm)
        end do
      end do
    end do

  case (2)
    ! in act_space
    if (jpad /= jpadl) return
    jph = mg1
    jpl = mg2
    !iwal = mg3
    !iwar = mg4
    iwupwei = jpad_upwei(jpad)
    isegdownwei = iseg_downwei(ipae)
    jpy = jphy(jph)
    in_ = ihy(jpy)

    call mma_allocate(lopu,4,loputmp,label='lopu')
    call jl_ne_jr(mp,jpl,jpr,mg3,mg4,lopu)
    do lp=1,mp
      iwl = lopu(1,lp)-1
      iwr = lopu(2,lp)-1
      jpe = lopu(3,lp)
      lwnu = iy(1,jpe)
      do jwu=jpy+1,jpy+in_
        iwal = iwl+ihy(jwu)
        iwar = iwr+ihy(jwu)
        do jwd=1,lwnu
          iwal = iwal+1
          iwar = iwar+1
          do iwd=0,iwupwei-1
            iwadl = iwalk_ad(jpadl,ipael,iwal,iwd)
            iwadr = iwalk_ad(jpad,ipae,iwar,iwd)
            do iwe=1,isegdownwei
              mm = iwadl+iwe
              nn = iwadr+iwe
              !vector2(mm) = vector2(mm)+vector1(nn)*wl
              !vector2(nn) = vector2(nn)+vector1(mm)*wl
              !vector2(mm) = vector2(mm)+vector1(nn)*wl
              !vector2(nn) = vector2(nn)+vector1(mm)*wl
              mg67 = ican_a(mg7)+mg6
              dm1tmp(mg67) = dm1tmp(mg67)+vector1(nn)*wl*vector1(mm)
              !write(nf2,'(4i4,4f18.10)') mg6,mg7,mm,nn,wl,vector1(mm),vector1(nn),dm1tmp(mg67)
            end do
          end do
        end do
      end do
    end do
    call mma_deallocate(lopu)

  case (3)
    ! between act and dbl
    jpl = mg1
    iwdl = mg2
    iwdr = mg3
    isegdownwei = iseg_downwei(ipae)

    call mma_allocate(lopu,4,loputmp,label='lopu')
    call jl_ne_jr(mp,jpl,jpr,mg4,mg5,lopu)
    do lp=1,mp
      iwal = lopu(1,lp)-1
      iwar = lopu(2,lp)-1
      jpe = lopu(3,lp)
      jwnu = iy(1,jpe)
      do ii=1,jwnu
        iwal = iwal+1
        iwar = iwar+1
        mm = iwalk_ad(jpadl,ipael,iwal,iwdl)
        nn = iwalk_ad(jpad,ipae,iwar,iwdr)
        do iwe=1,isegdownwei
          mm = mm+1        ! iwl=iwalk_ad
          nn = nn+1        ! iwl=iwalk_ad
          !vector2(mm) = vector2(mm)+vector1(nn)*wl
          !vector2(nn) = vector2(nn)+vector1(mm)*wl
          !vector2(mm) = vector2(mm)+vector1(nn)*wl
          !vector2(nn) = vector2(nn)+vector1(mm)*wl
          mg67 = ican_a(mg7)+mg6
          dm1tmp(mg67) = dm1tmp(mg67)+vector1(nn)*wl*vector1(mm)
          !write(nf2,'(a9,3i8,3f18.10)') '1_dbl_act',mg6,mm,nn,wl,vector1(mm),vector1(nn)
        end do
      end do
    end do
    call mma_deallocate(lopu)

end select

return

end subroutine prodab_h_1

subroutine prodab_h0_1(idb,mg1,mg2,mg3,mg4,mg5,wl,jpr,mg6,mg7)

use gugaci_global, only: dm1tmp, ican_a, ihy, ipae, ipael, iseg_downwei, iw_downwei, iy, jpad, jpad_upwei, jpadl, jphy, loputmp, &
                         nu_ae, vector1
use stdalloc, only: mma_allocate, mma_deallocate
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: idb, mg1, mg2, mg3, mg4, mg5, jpr, mg6, mg7
real(kind=wp), intent(in) :: wl
integer(kind=iwp) :: ii, in_, ipae_, isegdownwei, iwa, iwadl, iwadr, iwal, iwar, iwd, iwdl, iwdown, iwdr, iwe, iwl, iwr, iwupwei, &
                     jpe, jph, jpl, jpy, jwd, jwnu, jwu, lp, lwnu, m, mg67, mm, mp, nn
integer(kind=iwp), allocatable :: lopu(:,:)
integer(kind=iwp), external :: iwalk_ad

select case (idb)
  case default ! (1)
    ! in dbl_space
    jpad = mg2
    iwdl = mg3
    iwdr = mg4
    do ipae_=1,25
      ipae = ipae_ ! ipae is in global module, is this necessary?
      if (nu_ae(ipae) == 0) cycle
      iwdown = iw_downwei(jpad,ipae)
      if (iwdown == 0) cycle
      lwnu = iseg_downwei(ipae)
      do iwa=0,iwdown-1
        iwadl = iwalk_ad(jpad,ipae,iwa,iwdl)
        iwadr = iwalk_ad(jpad,ipae,iwa,iwdr)
        mm = iwadl
        nn = iwadr
        do m=1,lwnu
          mm = mm+1
          nn = nn+1
          !if (mm > nn) mntmp = mm*(mm-1)/2+nn
          !if (nn > mm) mntmp = nn*(nn-1)/2+mm
          !vector2(mntmp) = vector2(mntmp)+wl
          mg67 = ican_a(mg7)+mg6
          dm1tmp(mg67) = dm1tmp(mg67)+vector1(nn)*wl*vector1(mm)
          !if (mntmp == 2) write(u6,*) '  102',vector2(mntmp),wl
        end do
      end do
    end do

  case (2)
    ! in act_space
    if (jpad /= jpadl) return
    jph = mg1
    jpl = mg2
    !iwal = mg3
    !iwar = mg4
    iwupwei = jpad_upwei(jpad)
    isegdownwei = iseg_downwei(ipae)
    jpy = jphy(jph)
    in_ = ihy(jpy)

    call mma_allocate(lopu,4,loputmp,label='lopu')
    call jl_ne_jr(mp,jpl,jpr,mg3,mg4,lopu)
    do lp=1,mp
      iwl = lopu(1,lp)-1
      iwr = lopu(2,lp)-1
      jpe = lopu(3,lp)
      lwnu = iy(1,jpe)
      do jwu=jpy+1,jpy+in_
        iwal = iwl+ihy(jwu)
        iwar = iwr+ihy(jwu)
        do jwd=1,lwnu
          iwal = iwal+1
          iwar = iwar+1
          do iwd=0,iwupwei-1
            iwadl = iwalk_ad(jpadl,ipael,iwal,iwd)
            iwadr = iwalk_ad(jpad,ipae,iwar,iwd)
            do iwe=1,isegdownwei
              mm = iwadl+iwe
              nn = iwadr+iwe
              !if (mm > nn) mntmp = mm*(mm-1)/2+nn
              !if (nn > mm) mntmp = nn*(nn-1)/2+mm
              !vector2(mntmp) = vector2(mntmp)+wl
              mg67 = ican_a(mg7)+mg6
              dm1tmp(mg67) = dm1tmp(mg67)+vector1(nn)*wl*vector1(mm)
              !if (mntmp == 2) write(u6,*) '  202',vector2(mntmp),wl
            end do
          end do
        end do
      end do
    end do
    call mma_deallocate(lopu)

  case (3)
    ! between act and dbl
    jpl = mg1
    iwdl = mg2
    iwdr = mg3
    isegdownwei = iseg_downwei(ipae)

    call mma_allocate(lopu,4,loputmp,label='lopu')
    call jl_ne_jr(mp,jpl,jpr,mg4,mg5,lopu)
    do lp=1,mp
      iwal = lopu(1,lp)-1
      iwar = lopu(2,lp)-1
      jpe = lopu(3,lp)
      jwnu = iy(1,jpe)
      do ii=1,jwnu
        iwal = iwal+1
        iwar = iwar+1
        mm = iwalk_ad(jpadl,ipael,iwal,iwdl)
        nn = iwalk_ad(jpad,ipae,iwar,iwdr)
        do iwe=1,isegdownwei
          mm = mm+1             ! iwl=iwalk_ad
          nn = nn+1             ! iwl=iwalk_ad
          !if (mm > nn) mntmp = mm*(mm-1)/2+nn
          !if (nn > mm) mntmp = nn*(nn-1)/2+mm
          !vector2(mntmp) = vector2(mntmp)+wl
          mg67 = ican_a(mg7)+mg6
          dm1tmp(mg67) = dm1tmp(mg67)+vector1(nn)*wl*vector1(mm)
          !if (mntmp == 2) write(u6,*) '  302',vector2(mntmp),wl
        end do
      end do
    end do
    call mma_deallocate(lopu)

end select

return

end subroutine prodab_h0_1

!this subroutine prodab_2 does the dm2 part, which corresponds to vint_c
subroutine prodab_2(idb,mg1,mg2,mg3,mg4,mg5,wl,jpr,mg6)

use gugaci_global, only: log_prod
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: idb, mg1, mg2, mg3, mg4, mg5, jpr, mg6
real(kind=wp), intent(in) :: wl

if (log_prod == 1) call prodab_h_2(idb,mg1,mg2,mg3,mg4,mg5,wl,jpr,mg6)
if (log_prod == 2) call prodab_h0_2(idb,mg1,mg2,mg3,mg4,mg5,wl,jpr,mg6)

return

end subroutine prodab_2

subroutine prodab_h_2(idb,mg1,mg2,mg3,mg4,mg5,wl,jpr,mg6)

use gugaci_global, only: ihy, ipae, ipael, iseg_downwei, iw_downwei, iy, jpad, jpad_upwei, jpadl, jphy, loputmp, nu_ae, vector1, &
                         vector2
use stdalloc, only: mma_allocate, mma_deallocate
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: idb, mg1, mg2, mg3, mg4, mg5, jpr, mg6
real(kind=wp), intent(in) :: wl
integer(kind=iwp) :: ii, in_, ipae_, isegdownwei, iwa, iwadl, iwadr, iwal, iwar, iwd, iwdl, iwdown, iwdr, iwe, iwl, iwr, iwupwei, &
                     jpe, jph, jpl, jpy, jwd, jwnu, jwu, lp, lwnu, m, mm, mp, nn
integer(kind=iwp), allocatable :: lopu(:,:)
integer(kind=iwp), external :: iwalk_ad

select case (idb)
  case default ! (1)
    ! in dbl_space
    jpad = mg2
    iwdl = mg3
    iwdr = mg4
    do ipae_=1,25
      ipae = ipae_ ! ipae is in global module, is this necessary?
      if (nu_ae(ipae) == 0) cycle
      iwdown = iw_downwei(jpad,ipae)
      if (iwdown == 0) cycle
      lwnu = iseg_downwei(ipae)
      do iwa=0,iwdown-1
        iwadl = iwalk_ad(jpad,ipae,iwa,iwdl)
        iwadr = iwalk_ad(jpad,ipae,iwa,iwdr)
        mm = iwadl
        nn = iwadr
        do m=1,lwnu
          mm = mm+1
          nn = nn+1
          !vector2(mm) = vector2(mm)+vector1(nn)*wl
          !vector2(nn) = vector2(nn)+vector1(mm)*wl
          !vector2(mm) = vector2(mm)+vector1(nn)*wl
          !vector2(nn) = vector2(nn)+vector1(mm)*wl
          vector2(mg6) = vector2(mg6)+vector1(nn)*wl*vector1(mm)
          !if (mg6 == 29) write(nf2,'(i8,2i4,4f18.10)') mg6,mm,nn,vector2(mg6),vector1(mm),wl,vector1(nn)
        end do
      end do
    end do

  case (2)
    ! in act_space
    if (jpad /= jpadl) return
    jph = mg1
    jpl = mg2
    !iwal = mg3
    !iwar = mg4
    iwupwei = jpad_upwei(jpad)
    isegdownwei = iseg_downwei(ipae)
    jpy = jphy(jph)
    in_ = ihy(jpy)

    call mma_allocate(lopu,4,loputmp,label='loputmp')
    call jl_ne_jr(mp,jpl,jpr,mg3,mg4,lopu)
    do lp=1,mp
      iwl = lopu(1,lp)-1
      iwr = lopu(2,lp)-1
      jpe = lopu(3,lp)
      lwnu = iy(1,jpe)
      do jwu=jpy+1,jpy+in_
        iwal = iwl+ihy(jwu)
        iwar = iwr+ihy(jwu)
        do jwd=1,lwnu
          iwal = iwal+1
          iwar = iwar+1
          do iwd=0,iwupwei-1
            iwadl = iwalk_ad(jpadl,ipael,iwal,iwd)
            iwadr = iwalk_ad(jpad,ipae,iwar,iwd)
            do iwe=1,isegdownwei
              mm = iwadl+iwe
              nn = iwadr+iwe
              !vector2(mm) = vector2(mm)+vector1(nn)*wl
              !vector2(nn) = vector2(nn)+vector1(mm)*wl
              !vector2(mm) = vector2(mm)+vector1(nn)*wl
              !vector2(nn) = vector2(nn)+vector1(mm)*wl
              vector2(mg6) = vector2(mg6)+vector1(nn)*wl*vector1(mm)
              !if (mg6 == 15) write(nf2,'(2i4,4f18.10)') mm,nn,wl,vector1(mm),vector1(nn),vector2(mg6)
              !write(nf2,'(a3,2i4,i8,2f18.10)') 'act',mm,nn,mg6,vector1(mm),vector1(nn)
            end do
          end do
        end do
      end do
    end do
    call mma_deallocate(lopu)

  case (3)
    ! between act and dbl
    jpl = mg1
    iwdl = mg2
    iwdr = mg3
    isegdownwei = iseg_downwei(ipae)

    call mma_allocate(lopu,4,loputmp,label='loputmp')
    call jl_ne_jr(mp,jpl,jpr,mg4,mg5,lopu)
    do lp=1,mp
      iwal = lopu(1,lp)-1
      iwar = lopu(2,lp)-1
      jpe = lopu(3,lp)
      jwnu = iy(1,jpe)
      do ii=1,jwnu
        iwal = iwal+1
        iwar = iwar+1
        mm = iwalk_ad(jpadl,ipael,iwal,iwdl)
        nn = iwalk_ad(jpad,ipae,iwar,iwdr)
        do iwe=1,isegdownwei
          mm = mm+1             ! iwl=iwalk_ad
          nn = nn+1             ! iwl=iwalk_ad
          !vector2(mm) = vector2(mm)+vector1(nn)*wl
          !vector2(nn) = vector2(nn)+vector1(mm)*wl
          !vector2(mm) = vector2(mm)+vector1(nn)*wl
          !vector2(nn) = vector2(nn)+vector1(mm)*wl
          vector2(mg6) = vector2(mg6)+vector1(nn)*wl*vector1(mm)
        end do
      end do
    end do
    call mma_deallocate(lopu)

end select

return

end subroutine prodab_h_2

subroutine prodab_h0_2(idb,mg1,mg2,mg3,mg4,mg5,wl,jpr,mg6)

use gugaci_global, only: ihy, ipae, ipael, iseg_downwei, iw_downwei, iy, jpad, jpad_upwei, jpadl, jphy, loputmp, nu_ae, vector1, &
                         vector2
use stdalloc, only: mma_allocate, mma_deallocate
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: idb, mg1, mg2, mg3, mg4, mg5, jpr, mg6
real(kind=wp), intent(in) :: wl
integer(kind=iwp) :: ii, in_, ipae_, isegdownwei, iwa, iwadl, iwadr, iwal, iwar, iwd, iwdl, iwdown, iwdr, iwe, iwl, iwr, iwupwei, &
                     jpe, jph, jpl, jpy, jwd, jwnu, jwu, lp, lwnu, m, mm, mntmp, mp, nn
integer(kind=iwp), allocatable :: lopu(:,:)
integer(kind=iwp), external :: iwalk_ad

select case (idb)
  case default ! (1)
    ! in dbl_space
    jpad = mg2
    iwdl = mg3
    iwdr = mg4
    mntmp = 0
    do ipae_=1,25
      ipae = ipae_ ! ipae is in global module, is this necessary?
      if (nu_ae(ipae) == 0) cycle
      iwdown = iw_downwei(jpad,ipae)
      if (iwdown == 0) cycle
      lwnu = iseg_downwei(ipae)
      do iwa=0,iwdown-1
        iwadl = iwalk_ad(jpad,ipae,iwa,iwdl)
        iwadr = iwalk_ad(jpad,ipae,iwa,iwdr)
        mm = iwadl
        nn = iwadr
        do m=1,lwnu
          mm = mm+1
          nn = nn+1
          if (mm > nn) then
            mntmp = mm*(mm-1)/2+nn
          else
            mntmp = nn*(nn-1)/2+mm
          end if
          !vector2(mntmp) = vector2(mntmp)+wl
          vector2(mg6) = vector2(mg6)+vector1(mntmp)*wl*vector1(mntmp)

          !if (mntmp == 2) write(u6,*) '  102',vector2(mntmp),wl
        end do
      end do
    end do

  case (2)
    ! in act_space
    if (jpad /= jpadl) return
    jph = mg1
    jpl = mg2
    !iwal = mg3
    !iwar = mg4
    iwupwei = jpad_upwei(jpad)
    isegdownwei = iseg_downwei(ipae)
    jpy = jphy(jph)
    in_ = ihy(jpy)

    call mma_allocate(lopu,4,loputmp,label='loputmp')
    call jl_ne_jr(mp,jpl,jpr,mg3,mg4,lopu)
    do lp=1,mp
      iwl = lopu(1,lp)-1
      iwr = lopu(2,lp)-1
      jpe = lopu(3,lp)
      lwnu = iy(1,jpe)
      do jwu=jpy+1,jpy+in_
        iwal = iwl+ihy(jwu)
        iwar = iwr+ihy(jwu)
        do jwd=1,lwnu
          iwal = iwal+1
          iwar = iwar+1
          do iwd=0,iwupwei-1
            iwadl = iwalk_ad(jpadl,ipael,iwal,iwd)
            iwadr = iwalk_ad(jpad,ipae,iwar,iwd)
            do iwe=1,isegdownwei
              mm = iwadl+iwe
              nn = iwadr+iwe
              if (mm > nn) then
                mntmp = mm*(mm-1)/2+nn
              else
                mntmp = nn*(nn-1)/2+mm
              end if
              !vector2(mntmp) = vector2(mntmp)+wl
              vector2(mg6) = vector2(mg6)+vector1(mntmp)*wl*vector1(mntmp)

              !if (mntmp == 2) write(u6,*) '  202',vector2(mntmp),wl
            end do
          end do
        end do
      end do
    end do
    call mma_deallocate(lopu)

  case (3)
    ! between act and dbl
    jpl = mg1
    iwdl = mg2
    iwdr = mg3
    isegdownwei = iseg_downwei(ipae)

    call mma_allocate(lopu,4,loputmp,label='loputmp')
    call jl_ne_jr(mp,jpl,jpr,mg4,mg5,lopu)
    do lp=1,mp
      iwal = lopu(1,lp)-1
      iwar = lopu(2,lp)-1
      jpe = lopu(3,lp)
      jwnu = iy(1,jpe)
      do ii=1,jwnu
        iwal = iwal+1
        iwar = iwar+1
        mm = iwalk_ad(jpadl,ipael,iwal,iwdl)
        nn = iwalk_ad(jpad,ipae,iwar,iwdr)
        do iwe=1,isegdownwei
          mm = mm+1             ! iwl=iwalk_ad
          nn = nn+1             ! iwl=iwalk_ad
          if (mm > nn) then
            mntmp = mm*(mm-1)/2+nn
          else
            mntmp = nn*(nn-1)/2+mm
          end if
          !vector2(mntmp) = vector2(mntmp)+wl
          vector2(mg6) = vector2(mg6)+vector1(mntmp)*wl*vector1(mntmp)

          !if (mntmp == 2) write(u6,*) '  302',vector2(mntmp),wl
        end do
      end do
    end do
    call mma_deallocate(lopu)

end select

return

end subroutine prodab_h0_2
