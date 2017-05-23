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
      Subroutine LDF_Init(DoPairs,Verbose,irc)
C
C     Thomas Bondo Pedersen, June 2010.
C
C     Purpose: initialization of the LDF code.
C
C     Return codes:
C        irc=-1: user input error
C        irc=0: success
C        irc=1: internal error
C        irc=2: general unclassified error
C
C     Initialization must be considered failed if irc is non-zero.
C
      Implicit None
      Logical DoPairs
      Logical Verbose
      Integer irc
#include "localdf.fh"
#include "localdf_print.fh"
#include "WrkSpc.fh"

      Character*8 SecNam
      Parameter (SecNam='LDF_Init')

      Integer nTask
      Parameter (nTask=4)

      Character*17 Task(nTask)

      Logical Timing

      Integer nSym, nShell_Valence, nShell_Auxiliary
      Integer ip_T, l_T, iTask

      Real*8 tC0, tW0, tC1, tW1

      Integer  LDF_nBasSh
      External LDF_nBasSh

      Integer i, j
      Integer iTime
      iTime(i,j)=ip_T-1+2*(j-1)+i

      ! Init return code
      irc=0

      ! Check for symmetry (symmetry not implemented)
      Call Get_iScalar('nSym',nSym)
      If (nSym.ne.1) Then
         Write(6,'(A,A)')
     &   SecNam,': Local DF not implemented with symmetry!'
         irc=-1
         Return
      End If

      ! Init timing
      Timing=iPrint.ge.Inf_DetailedTiming
      If (Timing) Then
         l_T=2*nTask
         Call GetMem('LDFINIT','Allo','Real',ip_T,l_T)
         Call Cho_dZero(Work(ip_T),l_T)
      Else
         l_T=0
         ip_T=0
      End If

      ! Get a fresh start with Seward.
      ! Get number of valence and auxiliary shells.
      If (Timing) Call CWTime(tC0,tW0)
      Call LDF_CleanSheet(nShell_Valence,nShell_Auxiliary)
      If (Timing) Then
         Call CWTime(tC1,tW1)
         iTask=1
         Task(iTask)='Seward Init......'
         Work(iTime(1,iTask))=tC1-tC0
         Work(iTime(2,iTask))=tW1-tW0
      End If

      ! Set global basis and shell info
      ! nShell_Valence and nShell_Auxiliary are saved in
      ! localdf_bas.fh in this routine.
      If (Timing) Call CWTime(tC0,tW0)
      Call LDF_SetSh(nShell_Valence,nShell_Auxiliary,Verbose,irc)
      If (irc.ne.0) Then
         Write(6,'(A,A,I8)')
     &   SecNam,': LDF_SetSh returned code',irc
         irc=1
         Return
      End If
      If (Timing) Then
         Call CWTime(tC1,tW1)
         iTask=2
         Task(iTask)='Shell Info.......'
         Work(iTime(1,iTask))=tC1-tC0
         Work(iTime(2,iTask))=tW1-tW0
      End If

      ! Initilize atom info: shells, aux shells, and coordinates on each
      ! atom.
      If (Timing) Call CWTime(tC0,tW0)
      Call LDF_SetAtomInfo(Verbose,irc)
      If (irc.ne.0) Then
         Write(6,'(A,A,I8)')
     &   SecNam,': LDF_SetAtomInfo returned code',irc
         irc=1
         Return
      End If
      If (Timing) Then
         Call CWTime(tC1,tW1)
         iTask=3
         Task(iTask)='Atom Info........'
         Work(iTime(1,iTask))=tC1-tC0
         Work(iTime(2,iTask))=tW1-tW0
      End If

      ! Initialize atom pair info
      If (Timing) Call CWTime(tC0,tW0)
      If (DoPairs) Then
         Call LDF_SetAtomPairInfo(UseUniqueAtomPairs,Verbose,irc)
         If (irc.ne.0) Then
            Write(6,'(A,A,I8)')
     &      SecNam,': LDF_SetAtomPairInfo returned code',irc
            irc=1
            Return
         End If
      End If
      If (Timing) Then
         Call CWTime(tC1,tW1)
         iTask=4
         Task(iTask)='Atom Pair Info...'
         Work(iTime(1,iTask))=tC1-tC0
         Work(iTime(2,iTask))=tW1-tW0
      End If

      ! Print timing
      If (Timing) Then
         Write(6,'(/,A)')
     &   'Detailed Timing of LDF Initialization (CPU,Wall in s):'
         Do iTask=1,nTask
            Write(6,'(A17,1X,F7.1,1X,F7.1)')
     &      Task(iTask),Work(iTime(1,iTask)),Work(iTime(2,iTask))
         End Do
         Call xFlush(6)
         Call GetMem('LDFINIT','Free','Real',ip_T,l_T)
      End If

      End
