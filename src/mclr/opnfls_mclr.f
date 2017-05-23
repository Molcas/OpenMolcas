************************************************************************
* This file is part of OpenMolcas.                                     *
*                                                                      *
* OpenMolcas is free software; you can redistribute it and/or modify   *
* it under the terms of the GNU Lesser General Public License, v. 2.1. *
* OpenMolcas is distributed in the hope that it will be useful, but it *
* is provided "as is" and without any express or implied warranties.   *
* For more details see the full text of the license in the file        *
* LICENSE or in <http://www.gnu.org/licenses/>.                        *
************************************************************************
      Subroutine OpnFls_MCLR(iPL)
************************************************************************
*                                                                      *
*     Open files.                                                      *
*                                                                      *
*----------------------------------------------------------------------*
*                                                                      *
*                                                                      *
*----------------------------------------------------------------------*
*                                                                      *
*     history: none                                                    *
*                                                                      *
************************************************************************
      Implicit Real*8 (a-h,o-z)
#include "Files_mclr.fh"

#include "Input.fh"
      Character*8 Label
      Logical FoundTwoEls,DoCholesky
*---------------------------------------------------------------------*
*     Start                                                           *
*---------------------------------------------------------------------*
*---  open the JOBIPH file -------------------------------------------*
*     Call DaName(LuJob,FnJob)
      Call DaName(LuTemp,FnTemp)
*---  open the ORDINT file -------------------------------------------*
      Call f_Inquire(FnTwo,FoundTwoEls)
      Call DecideOnDirect(.True.,FoundTwoEls,Direct,DoCholesky)
      If (Direct) Then
         Write (6,*) 'OpnFls: No direct option in MCLR'
         Call Abend()
      Else
         If (.NOT.DoCholesky) Then
            If (iPL.ge.2) Write(6,*) 'Ordinary integral handling'
            iRc=-1
            iOpt=0
            Call OpnOrd(iRc,iOpt,FnTwo,LuTwo)
            If ( iRc.ne.0 ) Then
               Write (6,*) 'OpnFls: Error opening ORDINT'
               Call Abend()
            End If
         End If
      End If
      iRc=-1
      iOpt=0
      Call f_Inquire(FnMCK,McKinley)
      Call f_Inquire(FnPT2,PT2)
      If (McKinley) Then
*        Write(*,*) 'Calculating response on perturbation from mckinley'
         Call OpnMck(iRc,iOpt,FnMck,LuMck)
         If ( iRc.ne.0 ) Then
            Write (6,*) 'OpnFls: Error opening MCKINT'
            Call Abend()
         End If
         iRc=-1
         idum=0
         iOpt=0
         Label='SYMOP'
         call RdMck(irc,iopt,Label,idum,chirr,idum)
         If ( iRc.ne.0 ) Then
            Write (6,*) 'OpnFls: Error reading MCKINT'
            Write (6,'(A,A)') 'Label=',Label
            Call Abend()
         End If
*
      Else If (PT2) Then
       If (iPL.ge.2)
     &     Write(6,*) 'Calculating lagrange multipliers for CASPT2'
       Call DaName(LuPT2,FnPT2)
      Else
       If (iPL.ge.2) Then
       Write(6,*)'No ',FnPT2,' or ' ,FNMCK, ', I hope that is OK'
       Write(6,*)'Seward mode is assumed, reading perturbation from ',
     &           FnOne
       End If
      End If
*----------------------------------------------------------------------*
*     Exit                                                             *
*----------------------------------------------------------------------*
      Return
      End
