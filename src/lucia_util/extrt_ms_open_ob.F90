!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!                                                                      *
! Copyright (C) 2001, Jeppe Olsen                                      *
!***********************************************************************

!#define _DEBUGPRINT_
subroutine EXTRT_MS_OPEN_OB(IDET_OC,IDET_MS,IDET_OPEN_MS,NEL)
! A determinant IDET_OC, IDET_MS is given. Extract spinprojections
! for open shells
!
! Jeppe Olsen, December 2001

use Definitions, only: iwp
#ifdef _DEBUGPRINT_
use Definitions, only: u6
#endif

#include "intent.fh"

implicit none
integer(kind=iwp), intent(in) :: NEL, IDET_OC(NEL), IDET_MS(NEL)
integer(kind=iwp), intent(_OUT_) :: IDET_OPEN_MS(*)
integer(kind=iwp) :: IEL, IOPEN

IEL = 1
IOPEN = 0
! Loop over electrons
do
  if (IEL < NEL) then
    if (IDET_OC(IEL) /= IDET_OC(IEL+1)) then
      ! Single occupied orbital so
      IOPEN = IOPEN+1
      IDET_OPEN_MS(IOPEN) = IDET_MS(IEL)
      IEL = IEL+1
    else
      IEL = IEL+2
    end if
  else
    ! Last electron was not identical to previous, so
    ! necessarily single occupied.
    IOPEN = IOPEN+1
    IDET_OPEN_MS(IOPEN) = IDET_MS(IEL)
    IEL = IEL+1
  end if
  if (IEL > NEL) exit
end do

#ifdef _DEBUGPRINT_
write(u6,*) ' Input det, occ and ms'
call IWRTMA(IDET_OC,1,NEL,1,NEL)
call IWRTMA(IDET_MS,1,NEL,1,NEL)
write(u6,*) ' Number of open orbitals = ',IOPEN
write(u6,*) ' Output det : ms of open orbitals'
call IWRTMA(IDET_OPEN_MS,1,IOPEN,1,IOPEN)
#endif

end subroutine EXTRT_MS_OPEN_OB
