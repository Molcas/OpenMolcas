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
! Copyright (C) 1990, Roland Lindh                                     *
!               1990, IBM                                              *
!***********************************************************************

subroutine Expnd_i(Array,n,m)
!***********************************************************************
!                                                                      *
! Object: to do an in place expansion of a triagularized matrix.       *
!                                                                      *
! Called from: Cntrct                                                  *
!                                                                      *
! Calling    : DCopy  (ESSL)                                           *
!                                                                      *
!     Author: Roland Lindh, IBM Almaden Research Center, San Jose, CA  *
!             May '90                                                  *
!***********************************************************************

use Index_Functions, only: nTri_Elem
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: n, m
real(kind=wp), intent(inout) :: Array(m,n*n)
integer(kind=iwp) :: i, ii, ij, j, ji, nij

nij = nTri_Elem(n)
do i=n,1,-1
  do j=n,i+1,-1
    ji = n*(i-1)+j
    ij = n*(j-1)+i
    if (nij /= ij) Array(:,ij) = Array(:,nij)
    if (nij /= ji) Array(:,ji) = Array(:,nij)
    nij = nij-1
  end do
  ii = n*(i-1)+i
  if (nij /= ii) Array(:,ii) = Array(:,nij)
  nij = nij-1
end do

return

end subroutine Expnd_i
