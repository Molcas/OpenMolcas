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
  use mcpdft_output,only:lf
  use Fock_util_global,only:docholesky

  implicit none
#include "rasdim.fh"
#include "rasscf.fh"
#include "general.fh"
#include "warnings.h"

  logical,intent(out) :: DSCF
  logical :: file_exists
  integer :: iOpt,return_code

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
  if(return_code /= _RC_ALL_IS_WELL_) then
    write(lf,*) "MC-PDFT tried to open a file (ONEINT) containing"
    write(lf,*) "one-electron integrals, but failed. Something is"
    write(lf,*) "wrong with the file. Most probably it is simply"
    write(lf,*) "missing: please check. It should have been created"
    write(lf,*) "by SEWARD. Perhaps it is in the wrong directory?"
    call abend()
  endif

!---  open the ORDINT file --------------------------------------------*
  call f_Inquire('ORDINT',file_exists)
  Call DecideOnDirect(.True.,file_exists,DSCF,DoCholesky)
  If(.not. DSCF .And. .Not. DoCholesky) then
    return_code = -1
    iOpt = 0
    Call OpnOrd(return_code,iOpt,'ORDINT',LUINTA)
    If(return_code /= _RC_ALL_IS_WELL_) Then
      Write(lf,*) 'MC-PDFT tried to open a file (ORDINT) containing'
      Write(lf,*) 'two-electron integrals, but failed. Something'
      Write(lf,*) 'is wrong with the file. Most probably it is'
      Write(lf,*) 'simply missing: Please check. It should have'
      Write(lf,*) 'been created by SEWARD. Perhaps it is in the'
      Write(lf,*) 'wrong directory?'
      Call Abend()
    EndIf
  Else
    call f_Inquire('RUNFILE',file_exists)
    If(.not. file_exists) Then
      Write(LF,*) 'MC-PDFT tried to open a file (RUNFILE) containing'
      Write(LF,*) 'data from previous program steps. Something'
      Write(LF,*) 'is wrong with the file. Most probably it is'
      Write(LF,*) 'simply missing: Please check. It should have'
      Write(LF,*) 'been created by SEWARD.'
      Call Abend()
    EndIf
  EndIf
!---  open the file carrying the transfromed two-electron integrals ---*
  Call DaName(LUINTM,'TRAINT')

  Return
End
