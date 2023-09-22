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
! Copyright (C) 1996-2006, Thorstein Thorsteinsson                     *
!               1996-2006, David L. Cooper                             *
!***********************************************************************
!***********************************************************************
!*                                                                     *
!*  PERMVB := Permute orbitals in V1 according to IPERM.               *
!*  PERMCI := Permute orbitals in V1 according to IPERM.               *
!*                                                                     *
!*  V1 is either full CI vector (NDET), or just                        *
!*  VB determinants (NDETVB).                                          *
!*                                                                     *
!***********************************************************************

subroutine permvb_cvb(v1,iperm)
! Permutes orbitals in V1 according to IPERM.

use Definitions, only: wp, iwp

implicit none
#include "main_cvb.fh"
real(kind=wp) :: v1(*)
integer(kind=iwp) :: iperm(norb)
#include "WrkSpc.fh"
integer(kind=iwp) :: ialg, iv1, k1, k10, k11, k12, k13, k14, k15, k16, k17, k2, k5, k6, k7, k8, k9
logical(kind=iwp) :: vb
integer(kind=iwp), external :: mavailr_cvb, mstacki_cvb, mstackr_cvb

vb = .true.
k1 = mstacki_cvb((norb+1)*(nalf+1))
k2 = mstacki_cvb((norb+1)*(nbet+1))
k5 = mstacki_cvb(norb+1)
k6 = mstacki_cvb(norb+1)
k7 = mstacki_cvb(norb+1)
k8 = mstacki_cvb(norb)
k9 = mstacki_cvb(norb)
k10 = mstacki_cvb(norb)
k11 = mstacki_cvb(norb)
k12 = mstacki_cvb(norb)
k13 = mstacki_cvb(nda)
k14 = mstackr_cvb(nda)
k15 = mstacki_cvb(ndb)
k16 = mstackr_cvb(ndb)
if (vb) then
  k17 = mstackr_cvb(ndetvb)
else if (mavailr_cvb() >= ndet) then
  ialg = 1
  k17 = mstackr_cvb(ndet)
else
  ialg = 2
  k17 = mstackr_cvb(nda)
end if
if (vb) then
  call permvb2_cvb(v1,iperm,vb,iwork(ll(11)),iwork(ll(12)),iwork(k1),iwork(k2),iwork(k5),iwork(k6),iwork(k7),iwork(k8),iwork(k9), &
                   iwork(k10),iwork(k11),iwork(k12),iwork(k13),work(k14),iwork(k15),work(k16),work(k17),ialg)
else
  iv1 = nint(v1(1))
  !DLC iwork(maddr_r2i_cvb(iaddr_ci(iv1))) --> work(iaddr_ci(iv1))
  call permvb2_cvb(work(iaddr_ci(iv1)),iperm,vb,iwork(ll(11)),iwork(ll(12)),iwork(k1),iwork(k2),iwork(k5),iwork(k6),iwork(k7), &
                   iwork(k8),iwork(k9),iwork(k10),iwork(k11),iwork(k12),iwork(k13),work(k14),iwork(k15),work(k16),work(k17),ialg)
end if
call mfreei_cvb(k1)

return

end subroutine permvb_cvb
