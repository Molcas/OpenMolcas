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

use Definitions, only: wp, iwp

implicit none
integer(kind=iwp) :: IDM, LI, NSP
real(kind=wp) :: ARR1(IDM,*), ARR2(IDM,LI,*)
integer(kind=iwp) :: I, IJ, J

IJ = 0
call ZEROMA(ARR2(1,1,1),1,IDM)
do I=2,LI
  do J=1,I-1
    IJ = IJ+1
    call DCOPY_(IDM,ARR1(1,IJ),1,ARR2(1,I,J),1)
    call DCOPY_(IDM,ARR1(1,IJ),1,ARR2(1,J,I),1)
  end do
  call ZEROMA(ARR2(1,I,I),1,IDM)
end do
if (NSP < 0) then
  do I=1,LI
    call VNEG_CHT3(ARR2(1,1,I),1,ARR2(1,1,I),1,IDM*I)
  end do
end if

return

end subroutine EXPA2_UHF
