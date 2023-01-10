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

subroutine EXPA2_UHF(ARR1,IDM,LI,NSP,ARR2)
! THIS SUBROUTINE EXPANDS THE SECOND INDEX OF A MATRIX ARR1

use Index_Functions, only: nTri_Elem
use Constants, only: Zero
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: IDM, LI, NSP
real(kind=wp), intent(in) :: ARR1(IDM,nTri_Elem(LI-1))
real(kind=wp), intent(out) :: ARR2(IDM,LI,LI)
integer(kind=iwp) :: I, IJ, J

IJ = 0
do I=1,LI
  do J=1,I-1
    IJ = IJ+1
    ARR2(:,I,J) = ARR1(:,IJ)
    ARR2(:,J,I) = ARR1(:,IJ)
  end do
  ARR2(:,I,I) = Zero
end do
if (NSP < 0) then
  do I=1,LI
    ARR2(:,1:I,I) = -ARR2(:,1:I,I)
  end do
end if

return

end subroutine EXPA2_UHF
