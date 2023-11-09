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

subroutine compl_cvb(a,nvec,n)
! Creates orthogonal complement.
! On entry : A is square (NxN) and contains NVEC vectors.
! On exit  : A is a full matrix, NVEC first vectors are untouched,
! remaining orthonormal vectors span the orthogonal complement.

use stdalloc, only: mma_allocate, mma_deallocate
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: nvec, n
real(kind=wp), intent(inout) :: a(n,n)
integer(kind=iwp) :: i, ierr, imx, j
real(kind=wp) :: cmx, dum(1)
real(kind=wp), allocatable :: awrk(:,:), bwrk(:,:), cwrk(:), dwrk(:,:)
real(kind=wp), external :: ddot_

call mma_allocate(awrk,n,nvec,label='awrk')
call mma_allocate(dwrk,n,n,label='dwrk')
awrk(:,:) = a(:,1:nvec)
call unitmat(dwrk,n)
call schmidt_cvb(awrk,nvec,dum,n,0)
call schmidtd_cvb(awrk,nvec,dwrk,n,dum,n,0)
call mma_deallocate(awrk)
! Sort N vectors in order of decreasing norms
call mma_allocate(bwrk,n,n,label='bwrk')
call mma_allocate(cwrk,n,label='cwrk')
do i=1,n
  cwrk(i) = ddot_(n,dwrk(:,i),1,dwrk(:,i),1)
end do
do j=1,n
  cmx = cwrk(1)
  imx = 1
  do i=2,n
    if (cwrk(i) > cmx) then
      cmx = cwrk(i)
      imx = i
    end if
  end do
  cwrk(imx) = -real(j,kind=wp)
  bwrk(:,j) = dwrk(:,imx)
end do
call mma_deallocate(dwrk)
call schmidt_cvb(bwrk,n,dum,n,0)
! Extract N-NVEC remaining vectors with largest norms
do i=1,n
  cwrk(i) = ddot_(n,bwrk(:,i),1,bwrk(:,i),1)
end do
do j=1,n-nvec
  cmx = cwrk(1)
  imx = 1
  do i=2,n
    if (cwrk(i) > cmx) then
      cmx = cwrk(i)
      imx = i
    end if
  end do
  cwrk(imx) = -real(j,kind=wp)
  a(:,nvec+j) = bwrk(:,imx)
end do
ierr = 0
call nize_cvb(a(:,nvec+1:),n-nvec,dum,n,0,ierr)
call mma_deallocate(bwrk)
call mma_deallocate(cwrk)

return

end subroutine compl_cvb
