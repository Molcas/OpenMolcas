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
!***********************************************************************

subroutine OpnFls_RASSCF(DSCF,DoCholesky)
!***********************************************************************
!                                                                      *
!     Open files.                                                      *
!                                                                      *
!----------------------------------------------------------------------*
!                                                                      *
!     written by:                                                      *
!     M.P. Fuelscher                                                   *
!     University of Lund, Sweden, 1993                                 *
!                                                                      *
!----------------------------------------------------------------------*
!                                                                      *
!     history: none                                                    *
!                                                                      *
!***********************************************************************

use general_data, only: JOBIPH, JOBOLD, LUONEL, LUINTA, LUINTM, LUQUNE, LUDAVID, ITERFILE
use Definitions, only: u6

implicit none
logical DSCF, DoCholesky
logical test
integer iOpt, iRC
integer, external :: IsFreeUnit

!----------------------------------------------------------------------*
!     Start                                                            *
!----------------------------------------------------------------------*
!---  define logical unit numbers -------------------------------------*
! Molecular orbital input file  this variable is not used
! File is opened and closed i.e. around calls to rdvec.
!LUStartOrb = 19
! Job interface unit (-1 shows it has not been opened!)
JOBIPH = -1
! Old RASSCF job interface for input MO's and CI vector
JOBOLD = -1
! AO one-electron integrals
LUONEL = 16
! AO two-electron integrals
LUINTA = 40
! MO two-electron integrals
LUINTM = 13
! Temporary unit used for QUNE update
LUQUNE = 27
! Temporary unit for diagonalization
LUDAVID = 37
! Opening the JOBIPH file is delayed till after input processing at end
! of READIN_RASSCF. Only then is file name known.

!---  open the ORDINT file --------------------------------------------*
call f_Inquire('ORDINT',test)
call DecideOnDirect(.true.,test,DSCF,DoCholesky)
if ((.not. DSCF) .and. (.not. DoCholesky)) then
  iRc = -1
  iOpt = 0
  call OpnOrd(iRc,iOpt,'ORDINT',LUINTA)
  if (iRc /= 0) then
    write(u6,*) 'RASSCF tried to open a file (ORDINT) containing'
    write(u6,*) 'two-electron integrals, but failed. Something'
    write(u6,*) 'is wrong with the file. Most probably it is'
    write(u6,*) 'simply missing: Please check. It should have'
    write(u6,*) 'been created by SEWARD. Perhaps it is in the'
    write(u6,*) 'wrong directory?'
    call Abend()
  end if
else
  call f_Inquire('RUNFILE',test)
  if (.not. test) then
    write(u6,*) 'RASSCF tried to open a file (RUNFILE) containing'
    write(u6,*) 'data from previous program steps. Something'
    write(u6,*) 'is wrong with the file. Most probably it is'
    write(u6,*) 'simply missing: Please check. It should have'
    write(u6,*) 'been created by SEWARD.'
    call Abend()
  end if
end if
!---  open the file carrying the transfromed two-electron integrals ---*
call DaName(LUINTM,'TRAINT')
!---  open the DAVID file carrying temporary CI and sigma vectors -----*
!     Note the unit number is defined in the module general_data
call DaName(LuDavid,'TEMP01')
!---  open the file carrying the hessian update vectors ---------------*
call DaName(LuQune,'TEMP02')
!
! Open file for storage of information on CI-iterations
!
IterFile = IsFreeUnit(10)
call molcas_open(IterFile,'CIITER')
!----------------------------------------------------------------------*
!     Exit                                                             *
!----------------------------------------------------------------------*

end subroutine OpnFls_RASSCF
