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

integer function CHO_RS2F(LAB,ISHLAB,ISYMAB,IRED)
!
! Purpose: return index in reduced set IRED (1,2,3) of
!          element LAB in shell pair ISHLAB (sym. ISYMAB).
!          If not included in this reduced set, 0 is returned.

use ChoSwp, only: nnBstRSh, iiBstRSh, IndRed

implicit real*8(a-h,o-z)
#include "cholesky.fh"
character*8 SECNAM
parameter(SECNAM='CHO_RS2F')
integer K, K2

CHO_RS2F = 0

K = IIBSTR(ISYMAB,IRED)+IIBSTRSH(ISYMAB,ISHLAB,IRED)
K2 = K+NNBSTRSH(ISYMAB,ISHLAB,IRED)
if (IRED == 1) then
  do while ((K < K2) .and. (CHO_RS2F == 0))
    K = K+1
    if (INDRED(K,1) == LAB) CHO_RS2F = K
  end do
else if ((IRED == 2) .or. (IRED == 3)) then
  do while ((K < K2) .and. (CHO_RS2F == 0))
    K = K+1
    if (INDRED(INDRED(K,IRED),1) == LAB) CHO_RS2F = K
  end do
else
  call CHO_QUIT('IRED error in '//SECNAM,104)
end if

end function CHO_RS2F
