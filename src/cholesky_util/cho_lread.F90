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

function CHO_LREAD(ISYM,LWRK)
!
! Purpose: return a reasonable scratch space dimension for reading
!          previous vectors using cho_getvec.

use Cholesky, only: CHO_IOVEC, InfVec, nnBstR, NumCho, NVecRS1
use Definitions, only: iwp

implicit none
integer(kind=iwp) :: CHO_LREAD
integer(kind=iwp), intent(in) :: ISYM, LWRK
integer(kind=iwp) :: IRED, JRED, JVEC, LEN1, LEN2, LEN3, LMIN
integer(kind=iwp), parameter :: MNVECRS1 = 5

if (CHO_IOVEC == 1) then
  if ((NVECRS1(ISYM) < 1) .and. (NUMCHO(ISYM) > 0)) then
    NVECRS1(ISYM) = 1
    JVEC = 1
    IRED = INFVEC(JVEC,2,ISYM)
    do while (JVEC < NUMCHO(ISYM))
      JVEC = JVEC+1
      JRED = INFVEC(JVEC,2,ISYM)
      if (JRED == IRED) then
        NVECRS1(ISYM) = NVECRS1(ISYM)+1
      else
        JVEC = NUMCHO(ISYM)
      end if
    end do
  end if
  LEN1 = LWRK/3-1
  LEN2 = max(NVECRS1(ISYM),MNVECRS1)*NNBSTR(ISYM,1)
  LEN3 = min(LEN1,LEN2)
  LMIN = 2*NNBSTR(ISYM,1)
  CHO_LREAD = max(LEN3,LMIN)+1
else if ((CHO_IOVEC == 2) .or. (CHO_IOVEC == 3) .or. (CHO_IOVEC == 4)) then
  LEN1 = LWRK/3-1
  LMIN = 2*NNBSTR(ISYM,1)
  CHO_LREAD = max(LEN1,LMIN)+1
else
  CHO_LREAD = 2*NNBSTR(ISYM,1)
end if

end function CHO_LREAD
