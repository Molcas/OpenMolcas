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

use Constants, only: Zero
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: NDIM, MBLOCK, IFIL
real(kind=wp), intent(in) :: A(NDIM)
integer(kind=iwp) :: I, IDUM(1), IMZERO, ISTOP, MMBLOCK, NBACK, NBLOCK, NLABEL, NTRANS, START
real(kind=wp) :: XNORM
real(kind=wp), external :: ddot_
integer(kind=iwp), parameter :: ICRAY = 1, IPACK = 1

!write(u6,*) ' entering TODSC'
if (IPACK /= 0) then
  ! Check norm of A before writing
  XNORM = ddot_(nDim,A,1,A,1)
  if (XNORM == Zero) then
    IMZERO = 1
  else
    IMZERO = 0
  end if
  !write(u6,*) ' I am going to call ITODS'
  MMBLOCK = MBLOCK
  if (MMBLOCK > 1) MMBLOCK = 1
  IDUM(1) = IMZERO
  call ITODS(IDUM,1,MMBLOCK,IFIL)
  !write(u6,*) ' back from ITODS'
  if (IMZERO == 1) return
end if

if ((MBLOCK >= 0) .or. (ICRAY == 1)) then

  NBLOCK = MBLOCK
  if (MBLOCK <= 0) NBLOCK = NDIM
  ISTOP = 0
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
    START = ISTOP+1
    ISTOP = START+NBLOCK-1
    NBACK = NBACK-NTRANS
    write(IFIL) (A(I),I=START,ISTOP),NLABEL
    if (NBACK == 0) exit
  end do
end if

if ((ICRAY == 0) .and. (MBLOCK < 0) .and. (NDIM > 0)) then
  !call SQFILE(IFIL,1,A,2*NDIM)
  call SysHalt('todsc')
end if

!write(u6,*) ' leaving TODSC'

return

end subroutine TODSC_MCLR
