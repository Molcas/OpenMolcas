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

subroutine ddrestart_cvb(c,axc,vec,hp,solp,maxdav,n,nvguess1,nvrestart1)

use casvb_global, only: ifollow, nroot
use stdalloc, only: mma_allocate, mma_deallocate
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: maxdav, n
real(kind=wp), intent(inout) :: c(n,maxdav), axc(n,maxdav)
real(kind=wp), intent(out) :: vec(n)
real(kind=wp), intent(in) :: hp(maxdav,maxdav), solp(maxdav)
integer(kind=iwp), intent(out) :: nvguess1, nvrestart1
integer(kind=iwp) :: ir, ir_use
real(kind=wp), allocatable :: eigval(:), eigvec(:,:)

call mma_allocate(eigval,maxdav,label='eigval')
call mma_allocate(eigvec,maxdav,maxdav,label='eigvec')
eigvec(:,:) = hp(:,:)
call mxdiag_cvb(eigvec,eigval,maxdav)
call mma_deallocate(eigval)

nvrestart1 = 0
nvguess1 = 0
call mxatb_cvb(c,solp,n,maxdav,1,vec)
if (ifollow <= 2) then
  ! (Put lower-lying solutions in AxC :)
  do ir=1,nroot-1
    if (ifollow == 1) then
      ir_use = maxdav-ir+1
    else
      ir_use = ir
    end if
    call mxatb_cvb(c,eigvec(:,ir_use),n,maxdav,1,axc(:,ir+1))
  end do
  c(:,2:nroot) = axc(:,2:nroot)
end if
c(:,1) = vec(:)
call mma_deallocate(eigvec)

return

end subroutine ddrestart_cvb
