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
************************************************************************
      Subroutine OpnFls_RASSCF_m(DSCF,DoCholesky)
************************************************************************
*                                                                      *
*     Open files.                                                      *
*                                                                      *
*----------------------------------------------------------------------*
*                                                                      *
*     written by:                                                      *
*     M.P. Fuelscher                                                   *
*     University of Lund, Sweden, 1993                                 *
*                                                                      *
*----------------------------------------------------------------------*
*                                                                      *
*     history: none                                                    *
*                                                                      *
************************************************************************
      Implicit Real*8 (A-H,O-Z)
#include "rasdim.fh"
#include "rasscf.fh"
#include "general.fh"
#include "output_ras.fh"
      Parameter (ROUTINE='OPNFLS  ')
#include "davctl.fh"
#include "qnctl.fh"
      Logical DSCF,test,DoCholesky
*----------------------------------------------------------------------*
*     Start                                                            *
*----------------------------------------------------------------------*
      Call qEnter('OpnFls')
*---  define logical unit numbers -------------------------------------*
*...  Molecular orbital input file  this variable is not used
*...  File is opened and closed i.e. around calls to rdvec.
*      LUStartOrb=19
*...  Job interface unit (-1 shows it has not been opened!)
      JOBIPH=-1
*...  Old RASSCF job interface for input MO's and CI vector
      JOBOLD=-1
*...  AO one-electron integrals
      LUONEL=16
*...  AO two-electron integrals
      LUINTA=40
*...  MO two-electron integrals
      LUINTM=13
*...  Temporary unit used for QUNE update
      LUQUNE=27
*...  Temporary unit for diagonalization
      LUDAVID=37
*...  general purpose communication file COMFILE
*     Note: subr. GetInf uses unit 33 as logiacl unit
      LUCOM=30
* Opening the JOBIPH file is delayed till after input processing at end
* of READIN_RASSCF. Only then is file name known.

*---  open the ORDINT file --------------------------------------------*
      call f_Inquire('ORDINT',test)
      Call DecideOnDirect(.True.,test,DSCF,DoCholesky)
      If ( .not. DSCF .And. .Not.DoCholesky) then
        iRc=-1
        iOpt=0
        Call OpnOrd(iRc,iOpt,'ORDINT',LUINTA)
        If ( iRc.ne.0 ) Then
          Write(LF,*)'RASSCF tried to open a file (ORDINT) containing'
          Write(LF,*)'two-electron integrals, but failed. Something'
          Write(LF,*)'is wrong with the file. Most probably it is'
          Write(LF,*)'simply missing: Please check. It should have'
          Write(LF,*)'been created by SEWARD. Perhaps it is in the'
          Write(LF,*)'wrong directory?'
          Call Abend()
        End If
      Else
        call f_Inquire('RUNFILE',test)
        If ( .not.test ) Then
          Write(LF,*)'RASSCF tried to open a file (RUNFILE) containing'
          Write(LF,*)'data from previous program steps. Something'
          Write(LF,*)'is wrong with the file. Most probably it is'
          Write(LF,*)'simply missing: Please check. It should have'
          Write(LF,*)'been created by SEWARD.'
          Call Abend()
        End If
      End If
*---  open the file carrying the transfromed two-electron integrals ---*
      Call DaName(LUINTM,'TRAINT')
*---  open the DAVID file carrying temporary CI and sigma vectros -----*
*     Note the unit number is defined in the davctl.fh file
      Call DaName(LuDavid,'TEMP01')
*---  open the file carrying the hessian update vectors ---------------*
      Call DaName(LuQune,'TEMP02')
*
* Open file for storage of information on CI-iterations
*
      IterFile = IsFreeUnit(10)
      call molcas_open(IterFile,'CIITER')
c      Open(Unit=IterFile,Status='Unknown',Form='Formatted',
c     &     File='CIITER')
*----------------------------------------------------------------------*
*     Exit                                                             *
*----------------------------------------------------------------------*
      Call qExit('OpnFls')
      Return
      End
