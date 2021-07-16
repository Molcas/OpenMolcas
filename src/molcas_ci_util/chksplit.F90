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

use Definitions, only: iwp

implicit none
integer(kind=iwp) :: IERR
#include "rasdim.fh"
#include "output_ras.fh"
#include "general.fh"
#include "splitcas.fh"
#include "warnings.h"

!if (DoSplitCAS) then
IERR = 0
if (lRootSplit > nConf) IERR = 1
if (IERR == 1) then
  write(LF,*)
  write(LF,*) '******************** ERROR *********************'
  write(LF,'(1X,A)') 'Input Error:'
  write(LF,*) ' Root you are looking for is not reachable within'
  write(LF,*) ' the selected active space.'
  write(LF,*) ' Try to select a bigger active space!'
  write(LF,*) ' Root selected by user = ',lRootSplit
  write(LF,*) ' Root reachable = ',nConf
  write(LF,*) '************************************************'
  call Quit(_RC_INPUT_ERROR_)
end if

!if (NumSplit) then
IERR = 0
if (iDimBlockA < lRootSplit) IERR = 1
if (IERR == 1) then
  write(LF,*)
  write(LF,*) '******************** ERROR **********************'
  write(LF,*) 'Input Error: AA-Block selected is too small!'
  write(LF,'(1X,A,I5)') ' Root to be optimized :',lRootSplit
  write(LF,'(1X,A,I5)') ' AA-Block dimension   :',iDimBlockA
  write(LF,*) 'AA-Block must be always equal or greater than root.'
  write(LF,*) 'In a NUSP calculation increase iDimBlockA'
  write(LF,*) 'In a ENSP calculation increase the energy-gap'
  write(LF,*) 'In a PESP calculation increase the percentage'
  write(LF,*) '*************************************************'
  call Quit(_RC_INPUT_ERROR_)
  !iDimBlockA = lRootSplit
  !ircWar = 1
  !write(LF,'(1X,A,I8)') 'iDimBlockA has been reset to =',iDimBlockA
end if

IERR = 0
if (iDimBlockA > mxDimBlockA) IERR = 1
if (IERR == 1) then
  write(LF,*)
  write(LF,*) '***************** ERROR *****************'
  write(LF,'(1X,A,I6)') 'Input Error: Max dim. BlockA exceeded',mxDimBlockA
  write(LF,*) 'iDimBlockA selected by user = ',iDimBlockA
  write(LF,*) 'If you are running a NUSP calculation, please, decrease the value of iDimBlockA!'
  write(LF,*) 'If you are running a ENSP calculation, please, decrease the energy-gap!'
  write(LF,*) 'If you are running a PESP calculation, please, decrease the percentage!'
  write(LF,*) '************************************************'
  call Quit(_RC_INPUT_ERROR_)
end if
!end If
!^ End chk over NumSplit

IERR = 0
if (ThrSplit < min_ThrSplit) IERR = 1
if (IERR == 1) then
  write(LF,*)
  write(LF,*) '***************** ERROR *****************'
  write(LF,'(1X,A,I6)') 'Input Error: ThrSplit too small'
  write(LF,'(1X,A,I6)') 'minimum value ThrSplit = ',min_ThrSplit
  write(LF,*) 'ThrSplit selected by user = ',ThrSplit
  write(LF,*) '************************************************'
  call Quit(_RC_INPUT_ERROR_)
end if

!end if

return

end subroutine ChkSplit
