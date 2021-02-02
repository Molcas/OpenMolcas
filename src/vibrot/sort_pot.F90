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
! Copyright (C) 1981, Per-Olof Widmark                                 *
!***********************************************************************
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!                                                                      C
!  THIS ROUTINE SORTS THE VECTORS X AND Y, WHERE THE VECTOR X IS       C
!  USED AS THE KEYS.                                                   C
!                                                                      C
!  DATE: 81 04 14                                                      C
!                                                                      C
!  AUTHOR: P-O WIDMARK                                                 C
!                                                                      C
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

subroutine SORT_POT(X,Y,N)

use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: N
real(kind=wp), intent(inout) :: X(N), Y(N)
integer(kind=iwp) :: I, J, IND
real(kind=wp) :: SMALL, Z

do I=2,N
  IND = I-1
  SMALL = X(IND)
  do J=I,N
    if (X(J) <= SMALL) then
      IND = J
      SMALL = X(IND)
    end if
  end do
  Z = X(I-1)
  X(I-1) = X(IND)
  X(IND) = Z
  Z = Y(I-1)
  Y(I-1) = Y(IND)
  Y(IND) = Z
end do

return

end subroutine SORT_POT
