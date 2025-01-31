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

subroutine TODSCP(A,NDIM,MBLOCK,IFIL)
! TRANSFER ARRAY REAL*8  A(LENGTH NDIM) TO DISCFIL IFIL IN
! RECORDS WITH LENGTH NBLOCK.
!
! Packed version : Store only nonzero elements
! Small elements should be xeroed outside

use Constants, only: Zero
use lucia_data, only: IDISK

implicit none
integer NDIM, MBLOCK, IFIL
real*8 A(*)
!-jwk-cleanup integer START, STOP
real*8, external :: INPROD
integer ISCR(2), IDUMMY(1)
integer, parameter :: LPBLK = 50000
integer IPAK(LPBLK)
real*8 XPAK(LPBLK)
integer IPACK, IMZERO, MMBLOCK, IELMNT, LBATCH
real*8 XNORM

!write(6,*) ' entering TODSCP, file = ',IFIL
!call XFLUSH(6)
IPACK = 1
if (IPACK /= 0) then
  ! Check norm of A before writing
  XNORM = INPROD(A,A,NDIM)
  if (XNORM == Zero) then
    IMZERO = 1
  else
    IMZERO = 0
  end if
  MMBLOCK = MBLOCK
  if (MMBLOCK > 2) MMBLOCK = 2

  ISCR(1) = IMZERO
  ! Packing
  ISCR(2) = 1
  !call ITODS(ISCR,2,MMBLOCK,IFIL)
  call ITODS(ISCR,2,2,IFIL)
  if (IMZERO == 1) goto 1001
end if

! Loop over packed records of dimension LPBLK
IELMNT = 0
1000 continue
! The next LPBLK elements
LBATCH = 0
! Obtain next batch of elemnts
999 continue
if (NDIM >= 1) then
  IELMNT = IELMNT+1
  if (A(IELMNT) /= ZERO) then
    LBATCH = LBATCH+1
    IPAK(LBATCH) = IELMNT
    XPAK(LBATCH) = A(IELMNT)
  end if
end if
if ((LBATCH == LPBLK) .or. (IELMNT == NDIM)) goto 998
goto 999
! Send to DISC
998 continue
IDUMMY(1) = LBATCH
call IDAFILE(IFIL,1,IDUMMY,1,IDISK(IFIL))
if (LBATCH > 0) then
  call IDAFILE(IFIL,1,IPAK,LBATCH,IDISK(IFIL))
  call DDAFILE(IFIL,1,XPAK,LBATCH,IDISK(IFIL))
end if
if (IELMNT == NDIM) then
  IDUMMY(1) = -1
  call IDAFILE(IFIL,1,IDUMMY,1,IDISK(IFIL))
else
  IDUMMY(1) = 0
  call IDAFILE(IFIL,1,IDUMMY,1,IDISK(IFIL))
  goto 1000
end if
! End of loop over records of truncated elements
1001 continue

end subroutine TODSCP
