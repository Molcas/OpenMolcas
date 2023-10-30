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

subroutine oneexc_cvb(cfrom,cto,vij,diag,iPvb)

use casvb_global, only: iform_ci, ndet, norb, projcas
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero
use Definitions, only: wp, iwp, u6

implicit none
real(kind=wp), intent(in) :: cfrom(0:ndet)
real(kind=wp), intent(inout) :: cto(0:ndet), vij(*)
logical(kind=iwp), intent(in) :: diag
integer(kind=iwp), intent(in) :: iPvb
integer(kind=iwp) :: icfrom, icto, idens, nvij
real(kind=wp), allocatable :: vij2(:)

idens = 0
icfrom = nint(cfrom(0))
icto = nint(cto(0))

if (iform_ci(icfrom) /= 0) then
  write(u6,*) ' Unsupported format in ONEEXC/ONEDENS :',iform_ci(icfrom)
  call abend_cvb()
else if (iform_ci(icto) /= 0) then
  write(u6,*) ' Unsupported format in ONEEXC/ONEDENS :',iform_ci(icto)
  call abend_cvb()
end if

call oneexc2_cvb(cfrom(1:),cto(1:),vij,diag,idens,iPvb)

! If projcas and iPvb=0 we asume the normal density/1-ex. is required:
if (projcas .and. (iPvb /= 0)) then
  if (diag) then
    nvij = norb*norb
  else
    nvij = norb*(norb-1)
  end if
  call mma_allocate(vij2,nvij,label='vij2')
  if (idens == 0) then
    vij2(:) = -vij(1:nvij)
  else
    vij2(:) = Zero
  end if
  call oneexc2_cvb(cfrom(1:),cto(1:),vij2,diag,idens,3-iPvb)
  if (idens == 1) vij(1:nvij) = vij(1:nvij)-vij2(:)
  call mma_deallocate(vij2)
end if

return

end subroutine oneexc_cvb
