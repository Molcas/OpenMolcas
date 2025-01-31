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

subroutine NXTIJ(I,J,NI,NJ,IJSM,NONEW)
! An ordered pair (I,J) is given,I<=NI,J<=NJ
!
! Find next pair, if IJSM /= 0, I >= J

NONEW = 0
100 continue
if (I < NI) then
  I = I+1
else
  if (J < NJ) then
    I = 1
    J = J+1
  else
    NONEW = 1
    goto 101
  end if
end if
if ((IJSM /= 0) .and. (I < J)) goto 100
101 continue

NTEST = 0
if (NTEST /= 0) write(6,*) ' next (i,j) pair ',I,J

end subroutine NXTIJ
