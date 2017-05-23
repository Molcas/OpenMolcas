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
* Copyright (C) 2010, Thomas Bondo Pedersen                            *
************************************************************************
      Subroutine LDF_Final(TermInts,irc)
C
C     Thomas Bondo Pedersen, June 2010.
C
C     Finalize Local Density Fitting.
C
C     Return codes:
C        irc=0: success
C        irc=1: internal error
C
      Implicit None
      Logical TermInts
      Integer irc
#include "localdf.fh"

      Character*9 SecNam
      Parameter (SecNam='LDF_Final')

      Integer nErr

      nErr=0

      If (LDF_Run_Mode.eq.LDF_RUN_INTERNAL) Then
         ! Write atom pair info to disk for later use
         Call LDF_WriteAtomPairInfo(irc)
         If (irc.ne.0) Then
            Write(6,'(//,A,A,I8)')
     &      SecNam,': LDF_WriteAtomPairInfo returned code',irc
            nErr=nErr+1
         End If
      Else If (LDF_Run_Mode.ne.LDF_RUN_EXTERNAL) Then
         Call WarningMessage(2,SecNam//' improper run mode!')
         Write(6,'(A,I9)') 'Run mode=',LDF_Run_Mode
         Call LDF_Quit(1)
      End If

      ! Unset atom to atom pair map (if set)
      Call LDF_UnsetA2AP()

      ! Unset atom pair info
      Call LDF_UnsetAtomPairInfo(irc)
      If (irc.ne.0) Then
         Write(6,'(//,A,A,I8)')
     &   SecNam,': LDF_UnsetAtomPairInfo returned code',irc
         nErr=nErr+1
      End If

      ! Unset atom info
      Call LDF_UnsetAtomInfo(irc)
      If (irc.ne.0) Then
         Write(6,'(//,A,A,I8)')
     &   SecNam,': LDF_UnsetAtomInfo returned code',irc
         nErr=nErr+1
      End If

      ! Unset shell info
      Call LDF_UnsetSh(irc)
      If (irc.ne.0) Then
         Write(6,'(//,A,A,I8)')
     &   SecNam,': LDF_UnsetSh returned code',irc
         nErr=nErr+1
      End If

      ! Terminate integral setup and free shell tables
      If (TermInts) Then
         Call Term_Ints(.False.,.True.)
         Call Free_iSD()
      End If

      ! Set status on runfile
      Call LDF_SetStatusOnRunFile(LDF_Set-1)

      ! Set return code
      If (nErr.ne.0) Then
         irc=1
      Else
         irc=0
      End If

      End
