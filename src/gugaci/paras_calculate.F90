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

subroutine paras_calculate()

use gugaci_global, only: ibsm_ext, iesm_ext, iwt_orb_ext, iwt_sm_s_ext, jb_sys, jroute_sys, ng_sm, nlsm_dbl, nlsm_ext, spin
                         !, n_electron, norb_all, norb_dbl, norb_dz, norb_ext
use Symmetry_Info, only: Mul
use Definitions, only: iwp

implicit none
#include "Molcas.fh"
integer(kind=iwp) :: iaend, iaorb, iasta, ibend, iborb, ibsta, icnttmp, ijsm, isma, ismb, isumtmp, iwt_sm_sab(mxSym), iwttmp, ni, &
                     nij, nj
!integer(kind=iwp), parameter :: inlptb_new(156) = [ &
!                                                  ! 1  2  3  4  5  6  7  8  9  10  11  12 13
!                                                  ! a^r=1
!                                                   -1, 0, 0, 0, 0, 0,-7, 0,-9,-10,  0,-12, 0, &
!                                                  ! a^l=2
!                                                    0, 0, 0, 0, 0,-6, 0,-8, 0,  0,-11,  0, 0, &
!                                                  ! b_r=3
!                                                    4, 0, 0, 0, 0, 0, 0, 0, 0,  0,  0,  0, 0, &
!                                                  ! b_l=4
!                                                    5, 0, 0, 0, 0, 0, 0, 0, 0,  0,  0,  0, 0, &
!                                                  ! b^r=5
!                                                    0, 7, 8,10,11, 0, 0, 0, 0,  0,  0,  0, 0, &
!                                                  ! b^l=6
!                                                    0, 0, 9, 0,12, 0, 0, 0, 0,  0,  0,  0, 9, &
!                                                  ! c^'=7
!                                                    1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 3, &
!                                                  ! c^"=8
!                                                    1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12,13, &
!                                                  ! d_l^r=9
!                                                    6, 0, 0, 0, 0, 0, 0, 0, 0,  0,  0,  0, 0, &
!                                                  ! d^r^r=10
!                                                    0,-2, 0,-4, 0, 0, 0, 0, 0,  0,  0,  0, 0, &
!                                                  ! d^r^l=11
!                                                    0, 0,-3, 0,-5, 0, 0, 0, 0,  0,  0,  0, 0, &
!                                                  ! c^'" =12
!                                                    1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 0], &
!***********************************************************************
!   ar      =1 (+a^r)     drr     =2 (+d^rr)   drl     =3 (+d^rl)
!   arbr    =4 (+d^rr)    arbl    =5 (+d^rl)   ard_l^r =6 (+a^l)
!   drrb^r  =7 (+a^r)     drlb^r  =8 (+a^l)    drlb^l  =9 (+a^r)
!   arbrb^r =10 (+a^r)    arblb^r =11 (+a^l)   arblb^l =12 (+a^r)
!   drl     =13 (*)
!***********************************************************************
!                      inlptb_new(169) = [ -1,  0,  4,  5,  0,  0,  1,  1,  6,  0,  0,  1, &
!                                           0,  0,  0,  0,  7,  0,  2,  2,  0, -2,  0,  2, &
!                                           0,  0,  0,  0,  8,  9,  3,  3,  0,  0, -3,  3, &
!                                           0,  0,  0,  0, 10,  0,  4,  4,  0, -4,  0,  4, &
!                                           0,  0,  0,  0, 11, 12,  5,  5,  0,  0, -5,  5, &
!                                           0, -6,  0,  0,  0,  0,  6,  6,  0,  0,  0,  6, &
!                                          -7,  0,  0,  0,  0,  0,  7,  7,  0,  0,  0,  7, &
!                                           0, -8,  0,  0,  0,  0,  8,  8,  0,  0,  0,  8, &
!                                          -9,  0,  0,  0,  0,  0,  9,  9,  0,  0,  0,  9, &
!                                         -10,  0,  0,  0,  0,  0, 10, 10,  0,  0,  0, 10, &
!                                           0,-11,  0,  0,  0,  0, 11, 11,  0,  0,  0, 11, &
!                                         -12,  0,  0,  0,  0,  0, 12, 12,  0,  0,  0, 12, &
!                                           0,  0,  0,  0,  0,  9,  3, 13,  0,  0,  0,  0], &
!                      indord_plptype(12) = [0,0,0,1,1,3,3,1,1,2,2,2]   !severe_new_error_1

!ja_sys = int(n_electron*Half-spin)-norb_dz
jb_sys = int(spin+spin)
!jc_sys = norb_all-ja_sys-jb_sys
!jsm_sys = ns_sm
!write(u6,'(a9,3(1x,i3))') 'ja,jb,jc ',ja_sys,jb_sys,jc_sys
!if (jb_sys == jb_sys/2*2) then    ????? 11.2
!  w0_drl2_44 = v_sqtwo
!  w0_drl2_1122 = v_onevsqtwo
!else
!  w0_drl2_44 = -v_sqtwo
!  w0_drl2_1122 = v_onevsqtwo
!end if

if (jb_sys == 0) jroute_sys = 1
if (jb_sys == 1) jroute_sys = 2
if (jb_sys > 1) jroute_sys = 3

!do i=1,ng_sm
!  iwt_ext(i,1) = 0
!  iwt_ext(i,2) = nlsm_ext(i)
!end do
!iwt_ext(1,1) = 1
!
!do i=1,ng_sm
!  ni = nlsm_ext(i)
!  do j=i,ng_sm
!    nj = nlsm_ext(j)
!    if (i == j) then
!      nij = ni*(ni-1)/2
!    else
!      nij = ni*nj
!    end if
!    ijsm = Mul(i,j)
!    iwt_ext(ijsm,3) = iwt_ext(ijsm,3)+nij
!    iwt_ext(ijsm,4) = iwt_ext(ijsm,4)+nij
!  end do
!end do
!iwt_ext(1,4) = iwt_ext(1,4)+norb_ext
!
!do i=1,ng_sm
!  iwt_dbl(i,1) = 0
!  iwt_dbl(i,2) = nlsm_dbl(Mul(i,jsm_sys))
!  if (jroute_sys > 1) then
!    iwt_dbl(i,3) = nlsm_dbl(Mul(i,jsm_sys))
!  end if
!end do
!iwt_dbl(ns_sm,1) = 1
!
!do i=1,ng_sm
!  ni = nlsm_dbl(i)
!  do j=i,ng_sm
!    nj = nlsm_dbl(j)
!    if (i == j) then
!      nij = ni*(ni-1)/2
!    else
!      nij = ni*nj
!    end if
!    ijsm = Mul(Mul(i,j),ns_sm)
!    iwt_dbl(ijsm,4) = iwt_dbl(ijsm,4)+nij
!    if (jroute_sys > 2) then
!      iwt_dbl(ijsm,5) = iwt_dbl(ijsm,5)+nij
!    end if
!    iwt_dbl(ijsm,6) = iwt_dbl(ijsm,6)+nij
!    if (jroute_sys > 1) then
!      iwt_dbl(ijsm,6) = iwt_dbl(ijsm,6)+nij
!    end if
!  end do
!end do
!iwt_dbl(ns_sm,6) = iwt_dbl(ns_sm,6)+norb_dbl

do ismb=1,ng_sm
  iwt_sm_sab(ismb) = 0
end do
do ismb=1,ng_sm
  ni = nlsm_ext(ismb)
  ibsta = ibsm_ext(ismb)
  ibend = iesm_ext(ismb)
  do isma=1,ismb
    nj = nlsm_ext(isma)
    if (isma == ismb) then
      nij = ni*(ni-1)/2
    else
      nij = ni*nj
    end if
    ijsm = Mul(isma,ismb)
    !iwt_sm_sab(ijsm) = iwt_sm_sab(ijsm)+nij
    ! iwt_sm_sab function as a tmp array
    !iwt_sm_ext(isma,ismb) = iwt_sm_sab(ijsm)

    iasta = ibsm_ext(isma)
    iaend = iesm_ext(isma)
    if (ismb == isma) then
      ibsta = iasta+1
    end if
    iwttmp = iwt_sm_sab(ijsm)
    do iborb=ibsta,ibend
      do iaorb=iasta,min(iaend,iborb-1)
        iwttmp = iwttmp+1
        !if ((iaorb <= 0) .or. (iaorb > max_extorb)) write(u6,*)
        !if ((iborb <= 0) .or. (iborb > max_extorb)) write(u6,*)
        iwt_orb_ext(iaorb,iborb) = iwttmp
      end do
    end do

    iwt_sm_sab(ijsm) = iwt_sm_sab(ijsm)+nij
  end do
end do
iwt_sm_s_ext = iwt_sm_sab(1)

isumtmp = 0
do ismb=1,ng_sm
  icnttmp = iwt_sm_sab(ismb)
  iwt_sm_sab(ismb) = isumtmp
  isumtmp = isumtmp+icnttmp
end do
! imap_revorbtoorb_ext()
!do ismb=1,ng_sm
!  ibsta = ibsm_ext(ismb)
!  ibend = iesm_ext(ismb)
!  do isma=1,ismb
!    iasta = ibsm_ext(isma)
!    iaend = iesm_ext(isma)
!    if (isma == ismb) iaend = iaend-1
!
!    ismab = Mul(isma,ismb)
!    iwttmp = iwt_sm_ext(isma,ismb)+iwt_sm_sab(ismab)
!    do iaorb=iasta,iaend
!      do iborb=max(iaorb+1,ibsta),ibend
!        iwttmp = iwttmp+1
!        iwt_revorb_ext(iaorb,iborb) = iwttmp
!        imap_revorbtoorb_ext(iwttmp) = iwt_orb_ext(iaorb,iborb)
!      end do
!    end do
!  end do
!end do

do ismb=1,ng_sm
  iwt_sm_sab(ismb) = 0
end do
do ismb=1,ng_sm
  ni = nlsm_dbl(ismb)
  !ibsta = ibsm_dbl(ismb)
  !ibend = iesm_dbl(ismb)
  do isma=1,ismb
    nj = nlsm_dbl(isma)
    if (isma == ismb) then
      nij = ni*(ni-1)/2
    else
      nij = ni*nj
    end if
    ijsm = Mul(isma,ismb)
    !iwt_sm_dbl(isma,ismb) = iwt_sm_sab(ijsm)

    !iasta = ibsm_dbl(isma)
    !iaend = iesm_dbl(isma)
    !if (ismb == isma) then
    !  ibsta = iasta+1
    !end if
    !iwttmp = iwt_sm_sab(ijsm)
    !do iborb=ibsta,ibend
    !  do iaorb=iasta,min(iaend,iborb-1)
    !    iwttmp = iwttmp+1
    !    iwt_orb_dbl(iaorb,iborb) = iwttmp
    !  end do
    !end do

    iwt_sm_sab(ijsm) = iwt_sm_sab(ijsm)+nij
  end do
end do

!isumtmp = 0
!do ismb=1,ng_sm
!  icnttmp = iwt_sm_sab(ismb)
!  iwtsmsabtmp(ismb) = isumtmp
!  isumtmp = isumtmp+icnttmp
!end do
! imap_revorbtoorb_dbl()
!do ismb=1,ng_sm
!  ibsta = ibsm_dbl(ismb)
!  ibend = iesm_dbl(ismb)
!  do isma=1,ismb
!    iasta = ibsm_dbl(isma)
!    iaend = iesm_dbl(isma)
!    if (isma == ismb) iaend = iaend-1
!    ismab = Mul(isma,ismb)
!    iwttmp = iwt_sm_dbl(isma,ismb)+iwtsmsabtmp(ismab)
!    do iaorb=iasta,iaend
!      do iborb=max(iaorb+1,ibsta),ibend
!        iwttmp = iwttmp+1
!        iwt_revorb_dbl(iaorb,iborb) = iwttmp
!        imap_revorbtoorb_dbl(iwttmp) = iwt_orb_dbl(iaorb,iborb)
!      end do
!    end do
!  end do
!end do

!do ism=1,ng_sm
!  jsm = Mul(ism,ns_sm)
!  if (ism < jsm) then
!    itmp = iwt_sm_sab(jsm)
!    iwt_sm_sab(jsm) = iwt_sm_sab(ism)
!    iwt_sm_sab(ism) = itmp
!  end if
!end do
!
!iwt_sm_s_dbl = iwt_sm_sab(ns_sm)
!if (jroute_sys > 2) then
!  iwt_sm_s_dbl = iwt_sm_s_dbl*2
!end if

end subroutine paras_calculate
