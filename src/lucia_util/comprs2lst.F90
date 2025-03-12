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
! Copyright (C) 1996, Jeppe Olsen                                      *
!***********************************************************************

subroutine COMPRS2LST(I1,XI1,N1,I2,XI2,N2,NKIN,NKOUT)
! Two lists of excitations/annihilations/creations are given.
! Compress to common nonvanishing entries
!
! Jeppe Olsen, November 1996

use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: N1, N2, NKIN
integer(kind=iwp), intent(inout) :: I1(NKIN,N1), I2(NKIN,N2)
real(kind=wp), intent(inout) :: XI1(NKIN,N1), XI2(NKIN,N2)
integer(kind=iwp), intent(out) :: NKOUT
integer(kind=iwp) :: I1ACT, I2ACT, K

NKOUT = 0
do K=1,NKIN
  I1ACT = 0
  if (any(I1(K,:) /= 0)) I1ACT = 1
  I2ACT = 0
  if (any(I2(K,:) /= 0)) I2ACT = 1
  if ((I1ACT == 1) .and. (I2ACT == 1)) then
    NKOUT = NKOUT+1
    if (NKOUT /= K) then
      I1(NKOUT,:) = I1(K,:)
      XI1(NKOUT,:) = XI1(K,:)
      I2(NKOUT,:) = I2(K,:)
      XI2(NKOUT,:) = XI2(K,:)
    end if
  end if
end do

end subroutine COMPRS2LST
