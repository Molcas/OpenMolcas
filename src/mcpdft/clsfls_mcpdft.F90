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
! Copyright (C) 1993, Markus P. Fuelscher                              *
!               2024, Matthew R. Hennefarth                            *
!***********************************************************************

subroutine close_files_mcpdft()

use Fock_util_global, only: docholesky
use general_data, only: jobiph, jobold, luintm
use Definitions, only: iwp, u6

implicit none
#include "warnings.h"
integer(kind=iwp) :: iOpt, return_code

!---  close the JOBOLD file -------------------------------------------*
if ((JOBOLD > 0) .and. (JOBOLD /= JOBIPH)) then
  call DaClos(JOBOLD)
  JOBOLD = -1
else if (JOBOLD > 0) then
  JOBOLD = -1
end if
!---  close the JOBIPH file -------------------------------------------*
if (JOBIPH > 0) then
  call DaClos(JOBIPH)
  JOBIPH = -1
end if
!---  close the ORDINT file -------------------------------------------*
if (.not. DoCholesky) then
  return_code = -1
  call ClsOrd(return_code)
  if (return_code /= _RC_ALL_IS_WELL_) call WarningMessage(1,'Failed to close the ORDINT file.')
end if
!---  close the file carrying the transformed two-electron integrals --*
call DaClos(LUINTM)

!--- close the one-electorn integral file
return_code = -1
iOpt = 0
call clsone(return_code,iOpt)
if (return_code /= _RC_ALL_IS_WELL_) then
  write(u6,*) 'Error when trying to close the one-electron integral file.'
  call abend()
end if

end subroutine close_files_mcpdft
