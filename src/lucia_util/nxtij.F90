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

!#define _DEBUGPRINT_
subroutine NXTIJ(I,J,NI,NJ,IJSM,NONEW)
! An ordered pair (I,J) is given,I<=NI,J<=NJ
!
! Find next pair, if IJSM /= 0, I >= J

use Definitions, only: iwp
#ifdef _DEBUGPRINT_
use Definitions, only: u6
#endif

implicit none
integer(kind=iwp), intent(inout) :: I, J
integer(kind=iwp), intent(in) :: NI, NJ, IJSM
integer(kind=iwp), intent(out) :: NONEW

NONEW = 0
do
  if (I < NI) then
    I = I+1
  else
    if (J < NJ) then
      I = 1
      J = J+1
    else
      NONEW = 1
      exit
    end if
  end if
  if ((IJSM == 0) .or. (I >= J)) exit
end do

#ifdef _DEBUGPRINT_
write(u6,*) ' next (i,j) pair ',I,J
#endif

end subroutine NXTIJ
