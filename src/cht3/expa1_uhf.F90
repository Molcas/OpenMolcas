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

implicit none
real*8 ARR1, ARR2
integer IDM, LI, NSP, IJ, I, J, K
dimension ARR1(*), ARR2(LI,LI,*)

if (NSP > 0) then
  IJ = 1
  do K=1,IDM
    do I=1,LI
      call DCOPY_(I,ARR1(IJ),1,ARR2(I,1,K),LI)
      call DCOPY_(I,ARR1(IJ),1,ARR2(1,I,K),1)
      IJ = IJ+I
    end do
  end do
else
  IJ = 1
  do K=1,IDM
    ARR2(1,1,K) = 0.d0
    do I=2,LI
      ARR2(I,I,K) = 0.d0
      call DCOPY_(I-1,ARR1(IJ),1,ARR2(I,1,K),LI)
      do J=1,I-1
        ARR2(J,I,K) = -ARR1(IJ)
        IJ = IJ+1
      end do
    end do
  end do
end if

return

end subroutine EXPA1_UHF
