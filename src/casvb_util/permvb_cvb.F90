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

use casvb_global, only: iapr, ixapr, nda, ndet, ndetvb, norb
use stdalloc, only: mma_allocate, mma_deallocate
use Definitions, only: wp, iwp

implicit none
real(kind=wp), intent(inout) :: v1(0:ndet)
integer(kind=iwp), intent(in) :: iperm(norb)
integer(kind=iwp) :: ialg, mavailr, v2len
logical(kind=iwp) :: vb
real(kind=wp), allocatable :: v2(:)

vb = .true.
if (vb) then
  v2len = ndetvb
else
  call mma_maxDBLE(mavailr)
  if (mavailr >= ndet) then
    ialg = 1
    v2len = ndet
  else
    ialg = 2
    v2len = nda
  end if
end if
call mma_allocate(v2,v2len,label='v2')
if (vb) then
  call permvb2_cvb(v1,iperm,vb,iapr,ixapr,v2,ialg)
else
  call permvb2_cvb(v1(1:),iperm,vb,iapr,ixapr,v2,ialg)
end if

call mma_deallocate(v2)

return

end subroutine permvb_cvb
