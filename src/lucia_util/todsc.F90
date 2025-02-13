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

subroutine TODSC(A,NDIM,MBLOCK,IFIL)
! TRANSFER ARRAY REAL*8 A(LENGTH NDIM) TO DISCFIL IFIL IN
! RECORDS WITH LENGTH NBLOCK.

use lucia_data, only: IDISK
use Constants, only: Zero
use Definitions, only: wp, iwp

#include "intent.fh"

implicit none
integer(kind=iwp), intent(in) :: NDIM, MBLOCK, IFIL
real(kind=wp), intent(_IN_) :: A(NDIM)
integer(kind=iwp) :: IDUMMY(1), IMZERO, IPACK, ISCR(2), ISTOP, MMBLOCK, NBACK, NBLOCK, NLABEL, NTRANS, START
real(kind=wp) :: XNORM
real(kind=wp), external :: dDot_

IPACK = 1
if (IPACK /= 0) then
  ! Check norm of A before writing
  XNORM = dDot_(NDIM,A,1,A,1)
  if (XNORM == Zero) then
    IMZERO = 1
  else
    IMZERO = 0
  end if
  MMBLOCK = MBLOCK
  if (MMBLOCK > 2) MMBLOCK = 2

  ISCR(1) = IMZERO
  ! No packing
  ISCR(2) = 0
  call ITODS(ISCR,2,2,IFIL)
  if (IMZERO == 1) return
end if

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
  call DDAFILE(IFIL,1,A(START),ISTOP-START+1,IDISK(IFIL))
  IDUMMY(1) = NLABEL
  call IDAFILE(IFIL,1,IDUMMY,1,IDISK(IFIL))
  if (NBACK == 0) exit
end do

end subroutine TODSC
