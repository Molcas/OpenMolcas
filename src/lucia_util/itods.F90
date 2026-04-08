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

subroutine ITODS(IA,NDIM,MBLOCK,IFIL)
! TRANSFER ARRAY INTEGER IA(LENGTH NDIM) TO DISCFIL IFIL IN
! RECORDS WITH LENGTH NBLOCK.

use lucia_data, only: IDISK
use Definitions, only: iwp

#include "intent.fh"

implicit none
integer(kind=iwp), intent(in) :: NDIM, MBLOCK, IFIL
integer(kind=iwp), intent(_IN_) :: IA(NDIM)
integer(kind=iwp) :: IDUMMY(1), ISTOP, NBACK, NBLOCK, NLABEL, NTRANS, START

NBLOCK = MBLOCK

if (NBLOCK <= 0) NBLOCK = NDIM
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
  call IDAFILE(IFIL,1,IA(START),ISTOP-START+1,IDISK(IFIL))
  IDUMMY(1) = NLABEL
  call IDAFILE(IFIL,1,IDUMMY,1,IDISK(IFIL))
  if (NBACK == 0) exit
end do

end subroutine ITODS
