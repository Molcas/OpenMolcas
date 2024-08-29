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
! Copyright (C) 1991, Roland Lindh                                     *
!***********************************************************************

function lEmpty(Coeff,n2,ld2,m2)
!***********************************************************************
!                                                                      *
! Object: to set if partial or whole contraction matrix is empty.      *
!                                                                      *
! Called from: TwoEl                                                   *
!                                                                      *
! Calling    : None                                                    *
!                                                                      *
!     Author: Roland Lindh, Dept. of Theoretical Chemistry,            *
!             University of Lund, SWEDEN                               *
!             June '91                                                 *
!***********************************************************************

use Constants, only: Zero
use Definitions, only: wp, iwp

implicit none
logical(kind=iwp) :: lEmpty
integer(kind=iwp), intent(in) :: n2, ld2, m2
real(kind=wp), intent(in) :: Coeff(ld2,m2)
real(kind=wp) :: Temp
integer(kind=iwp) :: i, j

lEmpty = .true.
Temp = Zero
do i=1,n2
  do j=1,m2
    Temp = Temp+abs(Coeff(i,j))
  end do
end do

lEmpty = Temp == Zero

return

end function lEmpty
