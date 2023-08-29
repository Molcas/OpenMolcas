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

function CHO_RS2F(LAB,ISHLAB,ISYMAB,IRED)
!
! Purpose: return index in reduced set IRED (1,2,3) of
!          element LAB in shell pair ISHLAB (sym. ISYMAB).
!          If not included in this reduced set, 0 is returned.

use Cholesky, only: iiBstR, iiBstRSh, IndRed, nnBstRSh
use Definitions, only: iwp

implicit none
integer(kind=iwp) :: CHO_RS2F
integer(kind=iwp), intent(in) :: LAB, ISHLAB, ISYMAB, IRED
integer(kind=iwp) :: K, K2
character(len=*), parameter :: SECNAM = 'CHO_RS2F'

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
