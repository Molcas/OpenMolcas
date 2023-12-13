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

subroutine EXPA1_UHF(ARR1,IDM,LI,NSP,ARR2)
! THIS SUBROUTINE EXPANDS THE FIRST INDEX OF A MATRIX ARR1
!
! LI*(LI+1)/2

use Constants, only: Zero
use Definitions, only: wp, iwp

implicit none
real(kind=wp), intent(in) :: ARR1(*)
integer(kind=iwp), intent(in) :: IDM, LI, NSP
real(kind=wp), intent(out) :: ARR2(LI,LI,IDM)
integer(kind=iwp) :: I, IJ, K

if (NSP > 0) then
  IJ = 1
  do K=1,IDM
    do I=1,LI
      ARR2(I,1:I,K) = ARR1(IJ:IJ+I-1)
      ARR2(1:I,I,K) = ARR1(IJ:IJ+I-1)
      IJ = IJ+I
    end do
  end do
else
  IJ = 1
  do K=1,IDM
    do I=1,LI
      ARR2(I,1:I-1,K) = ARR1(IJ:IJ+I-2)
      ARR2(1:I-1,I,K) = -ARR1(IJ:IJ+I-2)
      ARR2(I,I,K) = Zero
      IJ = IJ+I-1
    end do
  end do
end if

return

end subroutine EXPA1_UHF
