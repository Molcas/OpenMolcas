************************************************************************
* This file is part of OpenMolcas.                                     *
*                                                                      *
* OpenMolcas is free software; you can redistribute it and/or modify   *
* it under the terms of the GNU Lesser General Public License, v. 2.1. *
* OpenMolcas is distributed in the hope that it will be useful, but it *
* is provided "as is" and without any express or implied warranties.   *
* For more details see the full text of the license in the file        *
* LICENSE or in <http://www.gnu.org/licenses/>.                        *
*                                                                      *
* Copyright (C) 1993, Markus P. Fuelscher                              *
*               1993, Per Ake Malmqvist                                *
************************************************************************
      Subroutine OpnFls_CASPT2
************************************************************************
C  purpose:
C  - initialize logical unit numbers
C  - open files
*----------------------------------------------------------------------*
C  written by:
C  M.P. Fuelscher and P. AA. Malmqvist
C  University of Lund, Sweden, 1993
************************************************************************
      Implicit Real*8 (A-H,O-Z)
      CHARACTER(2) CVEC,CMAT
*---------------------------------------------------------------------*
C  Start
*---------------------------------------------------------------------*
#include "rasdim.fh"
#include "caspt2.fh"
#include "output.fh"

      Call qEnter('OpnFls')
*---  define logical unit numbers ------------------------------------*
C  AO two-electron integrals
      LUINTA=20
C  AO one-electron integrals
      LUONEA=9
*
C  Used during solution of the caspt2 eqs
      LUSOLV=40
C  Used to hold S, B, and T matrices
      LUSBT =45
      CALL DANAME_MF_wa(LUSOLV,'LUSOLV')
      CALL DANAME_MF_wa(LUSBT ,'LUSBT ')
C  Half transformed integrals (uv|rs)
      LUHLF1=50
C  Half transformed integrals (uq|xs)
      LUHLF2=60
C  Half transformed integrals (uq|rt)
      LUHLF3=70
*
      CALL DANAME_MF_wa(LUHLF1,'LUHLF1')
      CALL DANAME_MF_wa(LUHLF2,'LUHLF2')
      CALL DANAME_MF_wa(LUHLF3,'LUHLF3')

* Disk-resident arrays for equation solving:
      LUDRA=30
      CALL DANAME_MF_wa(LUDRA,'DRARR')
      LUDRATOT=31
      CALL DANAME_MF_wa(LUDRATOT,'DRARRT')
*
C-SVC: assign logical units for RHS arrays and open files for writing
      DO IVEC=1,6
        LURHS(IVEC)=50+IVEC
        write(unit=CVEC, fmt='(I2.2)') IVEC
        CALL DANAME_MF_WA(LURHS(IVEC),'RHS_'//CVEC)
      END DO
C-SVC: assign logical units for SBT arrays and open files for writing
      DO IMAT=1,4
        LUH0T(IMAT)=60+IMAT
        write(unit=CMAT, fmt='(I2.2)') IMAT
        CALL DANAME_MF_WA(LUH0T(IMAT),'H0T_'//CMAT)
      END DO
C  Temporary unit with density matrices
      LUDMAT=90
      CALL DANAME_MF_wa(LUDMAT,'LUDMAT')

*---  open the files -------------------------------------------------*
C  Job interface
C      JOBIPH=15
C      CALL DANAME(JOBIPH,'JOBIPH')
C  A new JOBIPH file
C      JOBMIX=11
C      CALL DANAME(JOBMIX,'JOBMIX')
C  Temporary unit with excited CI expansions
      LUCIEX=10
      CALL DANAME_wa(LUCIEX,'LUCIEX')
C  Temporary unit with MO one-electron integrals
      LUONEM=16
      CALL DANAME_wa(LUONEM,'MOLONE')
C  Temporary unit with MO two-electron integrals (uv|xt)
      LUINTM=80
      CALL DANAME_MF_wa(LUINTM,'MOLINT')
C  AO one-electron integrals
      Call f_Inquire('ORDINT',Found2)
      Call DecideOnDirect(.False.,Found2,IfDirect,IfChol)
      If (IfChol) then
*        IF(IPRGLB.GE.USUAL) WRITE(6,*) 'This is a Cholesky CASPT2'
      else
*        IF(IPRGLB.GE.USUAL) WRITE(6,*) 'This is a conventional CASPT2'
        iRc=-1
        iOpt=0
        Call OpnOrd(iRc,iOpt,'ORDINT',LUINTA)
        If ( iRc.ne.0 ) Then
          WRITE(6,*)'OPNFLS Error: Failed to open the ORDINT file.'
          CALL ABEND()
        End If
      End If
*----------------------------------------------------------------------*
C  Exit
*----------------------------------------------------------------------*
      Call qExit('OpnFls')
      Return
      End
