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
      Subroutine LDF_ReplicateData(LuC_Local,irc)
C
C     Thomas Bondo Pedersen, July 2010.
C
C     Purpose: replicate all local DF data on all nodes in a real
C     parallel run. Do nothing in serial runs.
C
      Implicit None
      Integer LuC_Local
      Integer irc
#if defined (_MOLCAS_MPP_)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C PARALLEL CODE
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
#include "para_info.fh"
#include "WrkSpc.fh"
#include "ldf_atom_pair_info.fh"
#include "localdf_print.fh"

#if defined (_DEBUG_)
      Character*17 SecNam
      Parameter (SecNam='LDF_ReplicateData')
#endif

      Integer  LDF_DiskAddressOfC, LDF_AtomPair_DiagDim
      Integer  LDF_nBas_Atom, LDF_nBasAux_Pair
      External LDF_DiskAddressOfC, LDF_AtomPair_DiagDim
      External LDF_nBas_Atom, LDF_nBasAux_Pair

      Character*8 Label

      Integer n(1)

      Integer iAtomPair
      Integer ip, l
      Integer iAddr, jAddr
      Integer LuC
      Integer ip_C, l_C

      Logical Timing
      Real*8 tC0, tC1, tW0, tW1
      Real*8 tC2, tC3, tW2, tW3
      Real*8 tCIO, tWIO

      Integer i, j
      Integer ip_D
      Integer AP_Atoms, AP_1CLinDep, AP_2CFunctions
      Logical isUnique
      ip_D(i)=iWork(ip_AP_Diag-1+i)
      AP_Atoms(i,j)=iWork(ip_AP_Atoms-1+2*(j-1)+i)
      AP_1CLinDep(i,j)=iWork(ip_AP_1CLinDep-1+2*(j-1)+i)
      AP_2CFunctions(i,j)=iWork(ip_AP_2CFunctions-1+2*(j-1)+i)
      isUnique(i)=iWork(ip_AP_Unique-1+i).eq.i

      ! Init return code
      irc=0

      ! Only replicate data in truly parallel runs
      If (nProcs.gt.1 .and. Is_Real_Par()) Then
         ! Set timing flag
         Timing=iPrint.ge.Inf_DetailedTiming
C==========================
C        Replicate 1CLinDep
C==========================
         If (Timing) Call CWTime(tC0,tW0)
         Do iAtomPair=1,NumberOfAtomPairs
            If (isUnique(iAtomPair)) Then
               n(1)=AP_1CLinDep(1,iAtomPair)
               Call GAiGOp_SCAL(n(1),'max')
               If (n(1).gt.0) Then
                  l=3*n(1)
                  If (AP_1CLinDep(1,iAtomPair).gt.0) Then
#if defined (_DEBUG_)
                     If (LDF_DiskAddressOfC(iAtomPair).lt.0) Then
                        Call WarningMessage(1,
     &                            SecNam//': Parallelization error [1]')
                        irc=1
                        Return
                     End If
#endif
                     ip=AP_1CLinDep(2,iAtomPair)
#if defined (_DEBUG_)
                     If (ip.lt.1) Then
                        Call WarningMessage(1,
     &                            SecNam//': Parallelization error [2]')
                        irc=1
                        Return
                     End If
#endif
                  Else
                     Write(Label,'(A,I5.5)') '1CL',iAtomPair-1
                     Call GetMem(Label,'Allo','Inte',ip,l)
                     iWork(ip_AP_1CLinDep+2*(iAtomPair-1))=n(1)
                     iWork(ip_AP_1CLinDep+2*(iAtomPair-1)+1)=ip
                     Call iZero(iWork(ip),l)
                  End If
                  Call GAiGOp(iWork(ip),l,'+')
               End If
            End If
         End Do
         If (Timing) Then
            Call CWTime(tC1,tW1)
            Write(6,'(A,2F15.2,A)')
     &      'Time for 1CLinDep replication........',
     &      tC1-tC0,tW1-tW0,' seconds'
            tC0=tC1
            tW0=tW1
         End If
C=============================
C        Replicate 2CFunctions
C=============================
         Do iAtomPair=1,NumberOfAtomPairs
            If (isUnique(iAtomPair)) Then
               n(1)=AP_2CFunctions(1,iAtomPair)
               Call GAiGOp_SCAL(n(1),'max')
               If (n(1).gt.0) Then
                  l=4*n(1)
                  If (AP_2CFunctions(1,iAtomPair).gt.0) Then
#if defined (_DEBUG_)
                     If (LDF_DiskAddressOfC(iAtomPair).lt.0) Then
                        Call WarningMessage(1,
     &                            SecNam//': Parallelization error [3]')
                        irc=1
                        Return
                     End If
#endif
                     ip=AP_2CFunctions(2,iAtomPair)
#if defined (_DEBUG_)
                     If (ip.lt.1) Then
                        Call WarningMessage(1,
     &                            SecNam//': Parallelization error [4]')
                        irc=1
                        Return
                     End If
#endif
                  Else
                     Write(Label,'(A,I5.5)') '2CF',iAtomPair-1
                     Call GetMem(Label,'Allo','Inte',ip,l)
                     iWork(ip_AP_2CFunctions+2*(iAtomPair-1))=n(1)
                     iWork(ip_AP_2CFunctions+2*(iAtomPair-1)+1)=ip
                     Call iZero(iWork(ip),l)
                  End If
                  Call GAiGOp(iWork(ip),l,'+')
               End If
            End If
         End Do
         If (Timing) Then
            Call CWTime(tC1,tW1)
            Write(6,'(A,2F15.2,A)')
     &      'Time for 2CFunctions replication.....',
     &      tC1-tC0,tW1-tW0,' seconds'
            tC0=tC1
            tW0=tW1
         End If
C=========================================
C        Replicate updated diagonal blocks
C=========================================
         Do iAtomPair=1,NumberOfAtompairs
            If (isUnique(iAtomPair)) Then
               l=LDF_AtomPair_DiagDim(iAtomPair)
               ip=ip_D(iAtomPair)
               If (LDF_DiskAddressOfC(iAtomPair).lt.0) Then
                  Call Cho_dZero(Work(ip),l)
               End If
               Call GAdGOp(Work(ip),l,'+')
            End If
         End Do
         If (Timing) Then
            Call CWTime(tC1,tW1)
            Write(6,'(A,2F15.2,A)')
     &      'Time for diagonal replication........',
     &      tC1-tC0,tW1-tW0,' seconds'
            tC0=tC1
            tW0=tW1
         End If
C======================================
C        Replicate fitting coefficients
C======================================
         ! Open full coefficient file
         LuC=7
         Call DAName_MF_WA(LuC,'LDFC')
         ! Allocate memory for coefficients
         l_C=0
         Do iAtomPair=1,NumberOfAtomPairs
            If (isUnique(iAtomPair)) Then
               l_C=max(l_C,LDF_nBas_Atom(AP_Atoms(1,iAtomPair))
     &                    *LDF_nBas_Atom(AP_Atoms(2,iAtomPair))
     &                    *LDF_nBasAux_Pair(iAtomPair))
            End If
         End Do
         Call GetMem('RPL_C','Allo','Real',ip_C,l_C)
         ! Init disk address for coefficients
         iAddr=0
         ! Init I/O timing
         tCIO=0.0d0
         tWIO=0.0d0
         ! Replicate data one atom pair at a time
         Do iAtomPair=1,NumberOfAtomPairs
            If (isUnique(iAtomPair)) Then
               l=LDF_nBas_Atom(AP_Atoms(1,iAtomPair))
     &          *LDF_nBas_Atom(AP_Atoms(2,iAtomPair))
     &          *LDF_nBasAux_Pair(iAtomPair)
               jAddr=LDF_DiskAddressOfC(iAtomPair)
               If (jAddr.ge.0) Then
                  Call CWTime(tC2,tW2)
                  Call dDAFile(LuC_Local,2,Work(ip_C),l,jAddr)
                  Call CWTime(tC3,tW3)
                  tCIO=tCIO+(tC3-tC2)
                  tWIO=tWIO+(tW3-tW2)
               Else
                  Call Cho_dZero(Work(ip_C),l)
               End If
               Call GAdGOp(Work(ip_C),l,'+')
               iWork(ip_AP_DiskC-1+iAtomPair)=iAddr
               Call CWTime(tC2,tW2)
               Call dDAFile(LuC,1,Work(ip_C),l,iAddr)
               Call CWTime(tC3,tW3)
               tCIO=tCIO+(tC3-tC2)
               tWIO=tWIO+(tW3-tW2)
            End If
         End Do
         ! Deallocate coefficient memory
         Call GetMem('RPL_C','Free','Real',ip_C,l_C)
         ! Close full coefficient file
         Call DAClos(LuC)
         If (Timing) Then
            Call CWTime(tC1,tW1)
            Write(6,'(A,2F15.2,A)')
     &      'Time for coefficient replication.....',
     &      tC1-tC0,tW1-tW0,' seconds'
            Write(6,'(A,2F15.2,A)')
     &      '   (of which I/O required............',
     &      tCIO,tWIO,' seconds)'
            Call xFlush(6)
         End If
      End If
#else
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C SERIAL CODE
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      irc=0
c Avoid unused argument warnings
      If (.False.) Call Unused_integer(LuC_Local)
#endif
      End
