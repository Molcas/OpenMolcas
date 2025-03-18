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

subroutine FRMDSC_MCLR(ARRAY,NDIM,MBLOCK,IFILE,IMZERO)
! TRANSFER ARRAY FROM DISC FILE IFILE

implicit real*8(A-H,O-Z)
dimension ARRAY(*)
dimension IDUM(1)

IPACK = 1
if (IPACK /= 0) then
  ! Read if ARRAY is zero
  call IFRMDS(IDUM,1,MBLOCK,IFILE)
  IMZERO = IDUM(1)
  if (IMZERO == 1) then
    ARRAY(1:NDIM) = 0.0d0
    goto 1001
  end if
end if

ICRAY = 1

if ((MBLOCK >= 0) .or. (ICRAY == 1)) then
  NBLOCK = MBLOCK
  if (MBLOCK <= 0) NBLOCK = NDIM
  IREST = NDIM
  IBASE = 0
100 continue
  if (IREST > NBLOCK) then
    read(IFILE) (ARRAY(IBASE+I),I=1,NBLOCK)
    IBASE = IBASE+NBLOCK
    IREST = IREST-NBLOCK
  else
    read(IFILE) (ARRAY(IBASE+I),I=1,IREST)
    IREST = 0
  end if
  if (IREST > 0) goto 100
end if

if ((MBLOCK < 0) .and. (NDIM > 0) .and. (ICRAY == 0)) then
  !call SQFILE(IFILE,2,ARRAY,2*NDIM)
  call SysHalt('frmdsc')
end if

1001 continue

return

end subroutine FRMDSC_MCLR
