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

use stdalloc, only: mma_allocate, mma_deallocate
use Definitions, only: wp, iwp

implicit none
#include "main_cvb.fh"
real(kind=wp) :: v1(*)
integer(kind=iwp) :: iperm(norb)
#include "WrkSpc.fh"
integer(kind=iwp) :: ialg, iv1, v2len
logical(kind=iwp) :: vb
real(kind=wp), allocatable :: v2(:)
integer(kind=iwp), external :: mavailr_cvb

vb = .true.
if (vb) then
  v2len = ndetvb
else if (mavailr_cvb() >= ndet) then
  ialg = 1
  v2len = ndet
else
  ialg = 2
  v2len = nda
end if
call mma_allocate(v2,v2len,label='v2')
if (vb) then
  call permvb2_cvb(v1,iperm,vb,iwork(ll(11)),iwork(ll(12)),v2,ialg)
else
  iv1 = nint(v1(1))
  !DLC iwork(maddr_r2i_cvb(iaddr_ci(iv1))) --> work(iaddr_ci(iv1))
  call permvb2_cvb(work(iaddr_ci(iv1)),iperm,vb,iwork(ll(11)),iwork(ll(12)),v2,ialg)
end if

call mma_deallocate(v2)

return

end subroutine permvb_cvb
