!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!***********************************************************************

subroutine SQUARE(A,B,ICB,IRB,NROW)

use Index_Functions, only: nTri_Elem
use Definitions, only: wp, iwp

#include "intent.fh"

implicit none
integer(kind=iwp), intent(in) :: ICB, IRB, NROW
real(kind=wp), intent(in) :: A(*)
real(kind=wp), intent(_OUT_) :: B(*)
integer(kind=iwp) :: IC, ICOL, II, IND, IR, IROW

! PAM Sep 06: The two special cases here account for
! almost all calls of this code, and written such as
! not to store with non-zero stride--at least some
! rudimentary consideration for cache.
if (ICB == 1) then

  do IC=0,NROW-1
    II = nTri_Elem(IC)
    B(IC*IRB+1:IC*(IRB+1)+1) = A(II+1:II+IC+1)
  end do
  do IC=0,NROW-2
    do IR=IC+1,NROW-1
      B(1+IR+IC*IRB) = B(1+IC+IR*IRB)
    end do
  end do

else if (IRB == 1) then

  do IC=0,NROW-1
    II = nTri_Elem(IC)
    B(IC*ICB+1:IC*(ICB+1)+1) = A(II+1:II+IC+1)
  end do
  do IC=0,NROW-2
    do IR=IC+1,NROW-1
      B(1+IR+IC*ICB) = B(1+IC+IR*ICB)
    end do
  end do

else

  ! General and inefficient code:
  IND = 0
  do IROW=0,NROW-1
    do ICOL=0,IROW
      IND = IND+1
      B(1+IROW*ICB+ICOL*IRB) = A(IND)
      B(1+ICOL*ICB+IROW*IRB) = A(IND)
    end do
  end do

end if

return

end subroutine SQUARE
