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
! Copyright (C) 2010, Giovanni Li Manni                                *
!***********************************************************************

subroutine ChkSplit()
!***********************************************************************
!     SplitCAS Check for obvious errors or violation  of limits        *
!----------------------------------------------------------------------*
!     written by:                                                      *
!     Giovanni Li Manni (GLMJ)                                         *
!     University of Geneva, Switzerland, 2010                          *
!----------------------------------------------------------------------*
!     history: none                                                    *
!***********************************************************************

use general_data, only: nConf
use splitcas_data, only: iDimBlockA, lRootSplit, ThrSplit
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp) :: IERR
integer(kind=iwp), parameter :: MxDimBlockA = 2000
real(kind=wp), parameter :: min_ThrSplit = 1.0e-12_wp

#include "warnings.h"

!if (DoSplitCAS) then
IERR = 0
if (lRootSplit > nConf) IERR = 1
if (IERR == 1) then
  write(u6,*)
  write(u6,*) '******************** ERROR *********************'
  write(u6,'(1X,A)') 'Input Error:'
  write(u6,*) ' Root you are looking for is not reachable within'
  write(u6,*) ' the selected active space.'
  write(u6,*) ' Try to select a bigger active space!'
  write(u6,*) ' Root selected by user = ',lRootSplit
  write(u6,*) ' Root reachable = ',nConf
  write(u6,*) '************************************************'
  call Quit(_RC_INPUT_ERROR_)
end if

!if (NumSplit) then
IERR = 0
if (iDimBlockA < lRootSplit) IERR = 1
if (IERR == 1) then
  write(u6,*)
  write(u6,*) '******************** ERROR **********************'
  write(u6,*) 'Input Error: AA-Block selected is too small!'
  write(u6,'(1X,A,I5)') ' Root to be optimized :',lRootSplit
  write(u6,'(1X,A,I5)') ' AA-Block dimension   :',iDimBlockA
  write(u6,*) 'AA-Block must be always equal or greater than root.'
  write(u6,*) 'In a NUSP calculation increase iDimBlockA'
  write(u6,*) 'In a ENSP calculation increase the energy-gap'
  write(u6,*) 'In a PESP calculation increase the percentage'
  write(u6,*) '*************************************************'
  call Quit(_RC_INPUT_ERROR_)
  !iDimBlockA = lRootSplit
  !ircWar = 1
  !write(u6,'(1X,A,I8)') 'iDimBlockA has been reset to =',iDimBlockA
end if

IERR = 0
if (iDimBlockA > mxDimBlockA) IERR = 1
if (IERR == 1) then
  write(u6,*)
  write(u6,*) '***************** ERROR *****************'
  write(u6,'(1X,A,I6)') 'Input Error: Max dim. BlockA exceeded',mxDimBlockA
  write(u6,*) 'iDimBlockA selected by user = ',iDimBlockA
  write(u6,*) 'If you are running a NUSP calculation, please, decrease the value of iDimBlockA!'
  write(u6,*) 'If you are running a ENSP calculation, please, decrease the energy-gap!'
  write(u6,*) 'If you are running a PESP calculation, please, decrease the percentage!'
  write(u6,*) '************************************************'
  call Quit(_RC_INPUT_ERROR_)
end if
!end If
!^ End chk over NumSplit

IERR = 0
if (ThrSplit < min_ThrSplit) IERR = 1
if (IERR == 1) then
  write(u6,*)
  write(u6,*) '***************** ERROR *****************'
  write(u6,'(1X,A,I6)') 'Input Error: ThrSplit too small'
  write(u6,'(1X,A,I6)') 'minimum value ThrSplit = ',min_ThrSplit
  write(u6,*) 'ThrSplit selected by user = ',ThrSplit
  write(u6,*) '************************************************'
  call Quit(_RC_INPUT_ERROR_)
end if

!end if

return

end subroutine ChkSplit
