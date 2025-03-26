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

subroutine TODSC_MCLR(A,NDIM,MBLOCK,IFIL)
! TRANSFER ARRAY DOUBLE PRECISION  A(LENGTH NDIM) TO DISCFIL IFIL IN
! RECORDS WITH LENGTH NBLOCK.

implicit real*8(A-H,O-Z)
dimension A(1)
integer START, stop, IDUM(1)

!write(6,*) ' entering TODSC'
IPACK = 1
if (IPACK /= 0) then
  ! Check norm of A before writing
  XNORM = ddot_(nDim,A,1,A,1)
  if (XNORM == 0.0d0) then
    IMZERO = 1
  else
    IMZERO = 0
  end if
  !write(6,*) ' I am going to call ITODS'
  MMBLOCK = MBLOCK
  if (MMBLOCK > 1) MMBLOCK = 1
  IDUM(1) = IMZERO
  call ITODS(IDUM,1,MMBLOCK,IFIL)
  !write(6,*) ' back from ITODS'
  if (IMZERO == 1) return
end if

ICRAY = 1
if ((MBLOCK >= 0) .or. (ICRAY == 1)) then

  NBLOCK = MBLOCK
  if (MBLOCK <= 0) NBLOCK = NDIM
  stop = 0
  NBACK = NDIM
  ! LOOP OVER RECORDS
  do
    if (NBACK <= NBLOCK) then
      NTRANS = NBACK
      NLABEL = -NTRANS
    else
      NTRANS = NBLOCK
      NLABEL = NTRANS
    end if
    START = stop+1
    stop = START+NBLOCK-1
    NBACK = NBACK-NTRANS
    write(IFIL) (A(I),I=START,stop),NLABEL
    if (NBACK == 0) exit
  end do
end if

if ((ICRAY == 0) .and. (MBLOCK < 0) .and. (NDIM > 0)) then
  !call SQFILE(IFIL,1,A,2*NDIM)
  call SysHalt('todsc')
end if

!write(6,*) ' leaving TODSC'

return

end subroutine TODSC_MCLR
