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
!               1993, Per Ake Malmqvist                                *
!***********************************************************************

subroutine OpnFls_CASPT2()
!***********************************************************************
!  purpose:
!  - initialize logical unit numbers
!  - open files
!----------------------------------------------------------------------*
!  written by:
!  M.P. Fuelscher and P. AA. Malmqvist
!  University of Lund, Sweden, 1993
!***********************************************************************

use caspt2_global, only: LUCIEX, LUDMAT, LUDRA, LUDRATOT, LUH0T, LUHLF1, LUHLF2, LUHLF3, LUINTM, LUONEM, LURHS, LUSBT, LUSOLV
use caspt2_module, only: IfChol
use Definitions, only: iwp, u6

implicit none
integer(kind=iwp) :: iMat, iOpt, iRC, iVec, LUINTA
logical(kind=iwp) :: Found2, IfDirect
character(len=2) :: CMAT, CVEC

!----------------------------------------------------------------------*
!  Start
!----------------------------------------------------------------------*

!---  define logical unit numbers -------------------------------------*
! AO two-electron integrals
LUINTA = 20

! Used during solution of the caspt2 eqs
LUSOLV = 40
! Used to hold S, B, and T matrices
LUSBT = 45
call DANAME_MF_wa(LUSOLV,'LUSOLV')
call DANAME_MF_wa(LUSBT,'LUSBT ')
! Half transformed integrals (uv|rs)
LUHLF1 = 50
! Half transformed integrals (uq|xs)
LUHLF2 = 60
! Half transformed integrals (uq|rt)
LUHLF3 = 70

call DANAME_MF_wa(LUHLF1,'LUHLF1')
call DANAME_MF_wa(LUHLF2,'LUHLF2')
call DANAME_MF_wa(LUHLF3,'LUHLF3')

! Disk-resident arrays for equation solving:
LUDRA = 30
call DANAME_MF_wa(LUDRA,'DRARR')
LUDRATOT = 31
call DANAME_MF_wa(LUDRATOT,'DRARRT')

!-SVC: assign logical units for RHS arrays and open files for writing
do IVEC=1,8
  LURHS(IVEC) = 50+IVEC
  write(CVEC,'(I2.2)') IVEC
  call DANAME_MF_WA(LURHS(IVEC),'RHS_'//CVEC)
end do
!-SVC: assign logical units for SBT arrays and open files for writing
do IMAT=1,4
  LUH0T(IMAT) = 60+IMAT
  write(CMAT,'(I2.2)') IMAT
  call DANAME_MF_WA(LUH0T(IMAT),'H0T_'//CMAT)
end do
! Temporary unit with density matrices
LUDMAT = 90
call DANAME_MF_wa(LUDMAT,'LUDMAT')

!---  open the files --------------------------------------------------*
! Job interface
!JOBIPH = 15
!call DANAME(JOBIPH,'JOBIPH')
! A new JOBIPH file
!JOBMIX = 11
!call DANAME(JOBMIX,'JOBMIX')
! Temporary unit with excited CI expansions
LUCIEX = 10
call DANAME_wa(LUCIEX,'LUCIEX')
! Temporary unit with MO one-electron integrals
LUONEM = 16
call DANAME_wa(LUONEM,'MOLONE')
! Temporary unit with MO two-electron integrals (uv|xt)
LUINTM = 80
call DANAME_MF_wa(LUINTM,'MOLINT')
! AO one-electron integrals
call f_Inquire('ORDINT',Found2)
call DecideOnDirect(.false.,Found2,IfDirect,IfChol)
if (IfChol) then
  !if (IPRGLB >= USUAL) write(u6,*) 'This is a Cholesky CASPT2'
else
  !if (IPRGLB >= USUAL) writE(u6,*) 'This is a conventional CASPT2'
  iRc = -1
  iOpt = 0
  call OpnOrd(iRc,iOpt,'ORDINT',LUINTA)
  if (iRc /= 0) then
    write(u6,*) 'OPNFLS Error: Failed to open the ORDINT file.'
    call ABEND()
  end if
end if

!----------------------------------------------------------------------*
!  Exit
!----------------------------------------------------------------------*

end subroutine OpnFls_CASPT2
