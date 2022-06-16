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

subroutine Exp_2(Vector,n1,n2,Array,Fact)
!***********************************************************************
!                                                                      *
! Object: expand an array.                                             *
!                                                                      *
!     Author: Roland Lindh, Dept. of Theoretical Chemistry,            *
!             University of Lund, SWEDEN                               *
!             October '91                                              *
!***********************************************************************

use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: n1, n2
real(kind=wp), intent(out) :: Vector(n1,n2)
real(kind=wp), intent(in) :: Array(n2), Fact
integer(kind=iwp) :: i1

do i1=1,n1
  Vector(i1,:) = Array*Fact
end do

return

end subroutine Exp_2
