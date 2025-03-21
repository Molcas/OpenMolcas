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
      Subroutine ClsFls_MCLR()
************************************************************************
*                                                                      *
*     Open files.                                                      *
*                                                                      *
*----------------------------------------------------------------------*
*                                                                      *
*     written by:                                                      *
*     University of Lund, Sweden, 1992                                 *
*                                                                      *
*----------------------------------------------------------------------*
*                                                                      *
*     history: none                                                    *
*                                                                      *
************************************************************************
      use MCLR_Data, only: SA
      use MCLR_Data, only: FnMck,LuCSF2SD,LuJob,LuMck,LuQDat,LuTemp,
     &                      LuTri1
      use input_mclr, only: iMethod,TwoStep,RASSI
      Implicit None
      Logical DoCholesky
      Integer AixRm, iRC,iOpt
*---------------------------------------------------------------------*
*     Start                                                           *
*---------------------------------------------------------------------*
      If (iMethod.eq.2) Then
         Call DaClos(LuCSF2sd)
*------  close the JOBIPH file -------------------------------------------*
         Call DaClos(LuJob)
      End If
      Call DaClos(LuTemp)
*---  close the ORDINT file -------------------------------------------*
      Call DecideonCholesky(DoCholesky)
      If (.NOT.DoCholesky) then
         iRc=-1
         Call ClsOrd(iRc)
         If ( iRc.ne.0 ) Then
            Write (6,*) 'ClsFls: Error closing ORDINT'
            Call Abend()
         End If
      End If
      Call DaClos(LuTri1)
      If(TwoStep) Then
        Call DaClos(LuQDAT)
        !Call DaClos(LuMOTRA)
      End If
*
*---  Close the MckInt file or Remove the MCKINT file if SA---------------*
*     Do not remove file if we are producing data on the MckInt file for
*     the RASSI module!
*
      If (SA.and..Not.RASSI) Then
*        What the...? No control at all on what file is being removed!
*        call DaEras(LuMck)
         call DaClos(LuMck)
         iRC=AixRM(FnMck)
      Else
         iRc=-1
         iOpt=0
         Call ClsMck(iRc,iOpt)
         If ( iRc.ne.0 ) Then
            Write (6,*) 'ClsFls: Error closing MCKINT'
            Call Abend()
         End If
      End If
*
      Call ipTerm()
*----------------------------------------------------------------------*
*     Exit                                                             *
*----------------------------------------------------------------------*
      Return
      End
