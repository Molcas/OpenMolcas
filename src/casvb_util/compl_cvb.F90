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
integer(kind=iwp) :: nvec, n
real(kind=wp) :: a(n,n)
integer(kind=iwp) :: i, imx, j
real(kind=wp) :: cmx, dum(1)
real(kind=wp), allocatable :: awrk(:,:), bwrk(:,:), cwrk(:)
real(kind=wp), external :: ddot_

call mma_allocate(awrk,n,nvec+n,label='awrk')
call fmove_cvb(a,awrk,n*nvec)
call mxunit_cvb(awrk(1,1+nvec),n)
call schmidt_cvb(awrk,nvec,dum,n,0)
call schmidtd_cvb(awrk,nvec,awrk(1,nvec+1),n,dum,n,0)
! Sort N vectors in order of decreasing norms
call mma_allocate(bwrk,n,n,label='bwrk')
call mma_allocate(cwrk,n,label='cwrk')
do i=1,n
  cwrk(i) = ddot_(n,awrk(1,i+nvec),1,awrk(1,i+nvec),1)
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
  call fmove_cvb(awrk(1,imx+nvec),bwrk(1,j),n)
end do
call mma_deallocate(awrk)
call schmidt_cvb(bwrk,n,dum,n,0)
! Extract N-NVEC remaining vectors with largest norms
do i=1,n
  cwrk(i) = ddot_(n,bwrk(1,i),1,bwrk(1,i),1)
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
  call fmove_cvb(bwrk(1,imx),a(1,nvec+j),n)
end do
call nize_cvb(a(1,nvec+1),n-nvec,dum,n,0,0)
call mma_deallocate(bwrk)
call mma_deallocate(cwrk)

return

end subroutine compl_cvb
