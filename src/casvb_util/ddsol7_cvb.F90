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

subroutine ddsol7_cvb( &
#                     define _CALLING_
#                     include "ddsol_interface.fh"
                     )
! Solve linear equation in Davidson subspace.

use casvb_global, only: ifollow, ipdd, iroot, jroot, nfrdim1 => nfrdim, nroot
use stdalloc, only: mma_allocate, mma_deallocate
use Definitions, only: wp, iwp, u6

implicit none
#include "ddsol_interface.fh"
integer(kind=iwp) :: i
real(kind=wp) :: del, delmin
real(kind=wp), allocatable :: eigval(:), eigvec(:,:)

#include "macros.fh"
unused_var(rhsp)
unused_var(nfrdim)

if (ipdd >= 3) then
  write(u6,*) ' HP matrix (b) :'
  call mxprint2_cvb(ap,itdav,maxdav,itdav,0)
end if

call mma_allocate(eigval,itdav,label='eigval')
call mma_allocate(eigvec,itdav,itdav,label='eigvec')
eigvec(:,:) = ap(1:itdav,1:itdav)
call mxdiag_cvb(eigvec,eigval,itdav)

if (ifollow <= 2) then
  iroot = nroot
  jroot = mod(itdav,nroot)
  if (jroot == 0) jroot = nroot
  if ((itdav == maxdav) .or. (itdav == nfrdim1)) jroot = nroot
  iroot = min(itdav,iroot)
  jroot = min(itdav,jroot)
  if (ifollow == 1) then
    iroot = itdav-iroot+1
    jroot = itdav-jroot+1
  end if
else if (ifollow == 3) then
  write(u6,*) ' Overlap-based root following not yet implemented!'
  call abend_cvb()
else if (ifollow == 4) then
  ! Eigenvalue-based root following -- determine closest root:
  iroot = 1
  delmin = abs(eigval(1)-eig)
  do i=1,min(itdav,nroot)
    del = abs(eigval(i)-eig)
    if (del < delmin) then
      delmin = del
      iroot = i
    end if
  end do
  jroot = iroot
end if
eig = eigval(iroot)
solp(:) = eigvec(:,iroot)
eig_res = eigval(jroot)
solp_res(:) = eigvec(:,jroot)
if (ipdd >= 2) then
  write(u6,'(a)') ' Eigenvalues :'
  call vecprint_cvb(eigval,itdav)
  write(u6,'(a,i3,a)') ' Eigenvector number',iroot,' :'
  call vecprint_cvb(solp,itdav)
  if (jroot /= iroot) then
    write(u6,'(a,i3,a)') ' Eigenvector number',jroot,' :'
    call vecprint_cvb(solp_res,itdav)
  end if
end if
call mma_deallocate(eigval)
call mma_deallocate(eigvec)

return

end subroutine ddsol7_cvb
