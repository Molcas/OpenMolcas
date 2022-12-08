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
! Copyright (C) 2017, Roland Lindh                                     *
!***********************************************************************
! Version of Oct 21

subroutine XDIAXT(XDX,X,DIA,NDIM,SCR)
! Obtain XDX = X * DIA * X(Transposed)
! where DIA is a diagonal matrix stored as a vector

use Constants, only: Zero, One
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: NDIM
real(kind=wp), intent(out) :: XDX(NDIM,NDIM), SCR(NDIM,NDIM)
real(kind=wp), intent(in) :: X(NDIM,NDIM), DIA(NDIM)
integer(kind=iwp) :: I

! DIA * X
do I=1,NDIM
  SCR(:,I) = DIA(I)*X(:,I)
end do
! X * DIA * X(Transposed)
call dgemm_('N','T',NDIM,NDIM,NDIM,One,X,NDIM,SCR,NDIM,Zero,XDX,NDIM)

return

end subroutine XDIAXT
