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

subroutine span_cvb(a,nvec,nlin,s,n,metr)
! Creates span of a vector set.
! On entry : A is NVEC by N
! On exit  : NLIN is the number of linearly independent vectors
!            A contains set of NLIN linearly independent vectors
!              spanning the same set as the NVEC input vectors
!              Vectors will be orthonormal on exit

use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: nvec, n, metr
real(kind=wp), intent(inout) :: a(n,nvec)
integer(kind=iwp), intent(out) :: nlin
real(kind=wp), intent(in) :: s(*)
integer(kind=iwp) :: i, ierr, nvect
real(kind=wp) :: cnrm
real(kind=wp), parameter :: thresh = 1.0e-10_wp
real(kind=wp), external :: dnrm2_

nvect = nvec
ierr = 1
call nize_cvb(a,nvect,s,n,metr,ierr)
call schmidt_cvb(a,nvect,s,n,metr)
nlin = 0
do i=1,nvect
  cnrm = dnrm2_(n,a(:,i),1)
  if (cnrm > thresh) then
    nlin = nlin+1
    a(:,nlin) = a(:,i)
  end if
end do
ierr = 1
call nize_cvb(a,nlin,s,n,metr,ierr)

return

end subroutine span_cvb
