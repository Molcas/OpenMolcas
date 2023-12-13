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
! Copyright (C) 1989-1992,1998,1999, Roland Lindh                      *
!               1990, IBM                                              *
!               2000-2015, Valera Veryazov                             *
!***********************************************************************

subroutine Grid_it(iRun,ireturn)
!  iRun =1 normal run, 0=truncated from scf
!***********************************************************************
!                                                                      *
!  Object: Driver for evaluation MO values on a grid.                  *
!                                                                      *
!  Author: Roland Lindh, IBM Almaden Research Center, San Jose, CA     *
!          July '89 - May '90                                          *
!                                                                      *
!          Modified to gradient calculations September 1991 -          *
!          February 1992.                                              *
!                                                                      *
!          Modified to evaluating the spin density and spin density    *
!          gradients on a grid                                         *
!                                                                      *
!          Modified to interface with the MSI Cerius 2.                *
!          April 1998                                                  *
!                                                                      *
!          Modified at June- Sept 1999                                 *
!                                                                      *
!   This code was rewritten from scratch by V. Veryazov                *
!          Lund, 2000-2015                                             *
!***********************************************************************

use grid_it_globals, only: isUHF, levelprint, LuVal, LuVal_ab
use Definitions, only: iwp

implicit none
integer(kind=iwp), intent(in) :: iRun
integer(kind=iwp), intent(out) :: ireturn
#include "warnings.h"
integer(kind=iwp) :: nDiff
logical(kind=iwp) :: DoRys
character(len=1024) :: INPORB
!character(len=120) :: Lines(17)
integer(kind=iwp), external :: IPRINTLEVEL

! Prologue

levelprint = IPRINTLEVEL(-1)
if ((iRun == 0) .and. (levelprint < 3)) then
  levelprint = 0
  levelprint = IPRINTLEVEL(levelprint)
end if
if (iRun == 1) then
  call SetTim()

  !call bXML('GRID_IT')

  ! Get the work space size

end if

nDiff = 0
DoRys = .false.
call IniSew(DoRys,nDiff)

!---- Read the input

iReturn = 0
call Input_Grid_It(iRun,INPORB)
if (iReturn == _RC_INVOKED_OTHER_MODULE_) then
  ! take care to close files and release the potential memory...
  !close(LuOrb)
  close(LuVal)
  if (isUHF) close(LuVal_ab)
else

  ! Start computing the spin density and spin density gradient at the grid.

  call DrvMO(iRun,INPORB)

  !-----At the end of the calculation free all memory to check for
  !     corruption of the memory.

end if

call ClsSew()

! Epilogue

!if (iRun == 1) then
!  call eXML('GRID_IT')
!else
!  write(u6,*) 'Input file for molcasgv was generated'
!end if
!ireturn = 0

return

end subroutine Grid_it
