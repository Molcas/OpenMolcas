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
!               2024, Matthew R Hennefarth                             *
!***********************************************************************

subroutine open_files_mcpdft(DSCF)

use Fock_util_global, only: docholesky
use general_data, only: jobiph, jobold, luinta, luintm, luonel
use Definitions, only: iwp, u6

implicit none
logical(kind=iwp), intent(out) :: DSCF
#include "warnings.h"
integer(kind=iwp) :: iOpt, return_code
logical(kind=iwp) :: file_exists

!---  define logical unit numbers -------------------------------------*
!...  Job interface unit (-1 shows it has not been opened!)
JOBIPH = -1
!...  Old RASSCF job interface for input MO's and CI vector
JOBOLD = -1
!...  AO one-electron integrals
LUONEL = 16
!...  AO two-electron integrals
LUINTA = 40
!...  MO two-electron integrals
LUINTM = 13
! Opening the JOBIPH file is delayed till after input processing at end
! of READIN_RASSCF. Only then is file name known.

!--- open the one-electron integrals

return_code = -1
iOpt = 0
call OpnOne(return_code,iOpt,'ONEINT',LUONEL)
if (return_code /= _RC_ALL_IS_WELL_) then
  write(u6,*) 'MC-PDFT tried to open a file (ONEINT) containing'
  write(u6,*) 'one-electron integrals, but failed. Something is'
  write(u6,*) 'wrong with the file. Most probably it is simply'
  write(u6,*) 'missing: please check. It should have been created'
  write(u6,*) 'by SEWARD. Perhaps it is in the wrong directory?'
  call Abend()
end if

!---  open the ORDINT file --------------------------------------------*
call f_Inquire('ORDINT',file_exists)
call DecideOnDirect(.true.,file_exists,DSCF,DoCholesky)
if ((.not. DSCF) .and. (.not. DoCholesky)) then
  return_code = -1
  iOpt = 0
  call OpnOrd(return_code,iOpt,'ORDINT',LUINTA)
  if (return_code /= _RC_ALL_IS_WELL_) then
    write(u6,*) 'MC-PDFT tried to open a file (ORDINT) containing'
    write(u6,*) 'two-electron integrals, but failed. Something'
    write(u6,*) 'is wrong with the file. Most probably it is'
    write(u6,*) 'simply missing: Please check. It should have'
    write(u6,*) 'been created by SEWARD. Perhaps it is in the'
    write(u6,*) 'wrong directory?'
    call Abend()
  end if
else
  call f_Inquire('RUNFILE',file_exists)
  if (.not. file_exists) then
    write(u6,*) 'MC-PDFT tried to open a file (RUNFILE) containing'
    write(u6,*) 'data from previous program steps. Something'
    write(u6,*) 'is wrong with the file. Most probably it is'
    write(u6,*) 'simply missing: Please check. It should have'
    write(u6,*) 'been created by SEWARD.'
    call Abend()
  end if
end if
!---  open the file carrying the transfromed two-electron integrals ---*
call DaName(LUINTM,'TRAINT')

end subroutine open_files_mcpdft
