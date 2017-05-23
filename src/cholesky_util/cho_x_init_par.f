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
      SubRoutine Cho_X_Init_Par(irc,isDF)
C
C     Purpose: setup for parallel Cholesky/DF.
C
      Implicit None
      Integer irc
      Logical isDF

      If (isDF) Then
         Call Cho_X_Init_Par_DF(irc)
      Else
         Call Cho_X_Init_Par_GenBak()
         Call Cho_X_Init_Par_Cho(irc)
      End If

      End
      SubRoutine Cho_X_Init_Par_DF(irc)
C
C     Purpose: setup for parallel DF.
C
      Implicit None
      Integer irc

      Character*17 SecNam
      Parameter (SecNam = 'Cho_X_Init_Par_DF')

      Logical LocDbg
#if defined (_DEBUG_)
      Parameter (LocDbg = .True.)
#else
      Parameter (LocDbg = .False.)
#endif

#if defined (_MOLCAS_MPP_)
#include "para_info.fh"
#include "cholesky.fh"
      Integer nV(8)
      Integer iSym
      Logical isSerial

      irc = 0

C     Return if serial.
C     -----------------

      isSerial = nProcs.eq.1 .or. .not.Is_Real_Par()
      If (isSerial) Then
         If (LocDbg) Then
            Write(6,*) SecNam,': serial run, nothing to do...'
            Write(6,*) '#nodes: ',nProcs,'  myRank: ',myRank
         End If
         Return
      Else
         If (LocDbg) Then
            Write(6,*) SecNam,': parallel run...'
            Write(6,*) '#nodes: ',nProcs,'  myRank: ',myRank
         End If
      End If

C     Reset number of vectors to the number on this node as stored on
C     the runfile.
C     ---------------------------------------------------------------

      Call iCopy(nSym,NumCho,1,nV,1)
      Call Get_iArray('nVec_RI',NumCho,nSym)
      NumChT = NumCho(1)
      Do iSym = 2,nSym
         NumChT = NumChT + NumCho(iSym)
      End Do

C     Debug print.
C     ------------

      If (LocDbg) Then
         Write(6,*)
         Write(6,*) 'Output from ',SecNam,':'
         Write(6,*) 'NumCho before: ',(nV(iSym),iSym=1,nSym)
         Write(6,*) 'NumCho after : ',(NumCho(iSym),iSym=1,nSym)
      End If

#else

      irc = 0
      If (LocDbg) Then
         Write(6,*) SecNam,': serial run, nothing to do...'
      End If

#endif

      End
      SubRoutine Cho_X_Init_Par_Cho(irc)
C
C     Purpose: setup for parallel Cholesky.
C
      Implicit None
      Integer irc

      Character*18 SecNam
      Parameter (SecNam = 'Cho_X_Init_Par_Cho')

      Logical LocDbg
#if defined (_DEBUG_)
      Parameter (LocDbg = .True.)
#else
      Parameter (LocDbg = .False.)
#endif

#if defined (_MOLCAS_MPP_)
#include "para_info.fh"
#include "cholesky.fh"
#include "choptr.fh"
#include "WrkSpc.fh"

      Integer N2
      Parameter (N2 = InfVec_N2)

      Integer nV(8)
      Integer ip_IDV, l_IDV
      Integer ip_myInfV, l_myInfV
      Integer iSym
      Logical isSerial

      Integer InfVec, IDV, myInfV
      Integer i, j, k

      InfVec(i,j,k)=iWork(ip_InfVec-1+MaxVec*N2*(k-1)+MaxVec*(j-1)+i)
      IDV(i)=iWork(ip_IDV-1+i)
      myInfV(i)=iWork(ip_myInfV-1+i)

      irc = 0

C     Return if serial.
C     -----------------

      isSerial = nProcs.eq.1 .or. .not.Is_Real_Par()
      If (isSerial) Then
         If (LocDbg) Then
            Write(6,*) SecNam,': serial run, nothing to do...'
            Write(6,*) '#nodes: ',nProcs,'  myRank: ',myRank
         End If
         Return
      Else
         If (LocDbg) Then
            Write(6,*) SecNam,': parallel run...'
            Write(6,*) '#nodes: ',nProcs,'  myRank: ',myRank
         End If
      End If

C     Reset vector info to fit vectors stored on this node.
C     -----------------------------------------------------

      Do iSym = 1,nSym
         nV(iSym) = 0
         If (NumCho(iSym) .gt. 0) Then
            l_IDV = NumCho(iSym)
            Call GetMem('IDV','Allo','Inte',ip_IDV,l_IDV)
            Call Cho_Distrib_Vec(1,NumCho(iSym),iWork(ip_IDV),nV(iSym))
            If (nV(iSym) .gt. 0) Then
               l_myInfV = nV(iSym)
               Call GetMem('myInfV','Allo','Inte',ip_myInfV,l_myInfV)
               Do j = 1,N2
                  If (j .ne. 3) Then
                     k = ip_myInfV - 1
                     Do i = 1,nV(iSym)
                        iWork(k+i) = InfVec(IDV(i),j,iSym)
                     End Do
                     k = ip_InfVec - 1
     &                 + MaxVec*N2*(iSym-1) + MaxVec*(j-1)
                     Do i = 1,nV(iSym)
                        iWork(k+i) = myInfV(i)
                     End Do
                  End If
               End Do
               Call GetMem('myInfV','Free','Inte',ip_myInfV,l_myInfV)
            End If
            Call GetMem('IDV','Free','Inte',ip_IDV,l_IDV)
         End If
      End Do

C     Reset number of vectors.
C     ------------------------

      Call iSwap(nSym,NumCho,1,nV,1)
      NumChT = NumCho(1)
      Do iSym = 2,nSym
         NumChT = NumChT + NumCho(iSym)
      End Do

C     Debug print.
C     ------------

      If (LocDbg) Then
         Write(6,*)
         Write(6,*) 'Output from ',SecNam,':'
         Write(6,*) 'NumCho before: ',(nV(iSym),iSym=1,nSym)
         Write(6,*) 'NumCho after : ',(NumCho(iSym),iSym=1,nSym)
      End If

#else

      irc = 0
      If (LocDbg) Then
         Write(6,*) SecNam,': serial run, nothing to do...'
      End If

#endif

      End
      SubRoutine Cho_X_Init_Par_GenBak()
      Implicit None
#include "cholesky.fh"
#include "choptr.fh"
#include "chopar.fh"
#include "WrkSpc.fh"

#include "para_info.fh"

      If (Is_Real_Par()) Then
         l_InfVec_Bak = l_InfVec
         Call GetMem('InfVec_Bak','Allo','Inte',ip_InfVec_Bak,
     &                                           l_InfVec_Bak)
         Call iCopy(l_InfVec,iWork(ip_InfVec),1,iWork(ip_InfVec_Bak),1)
         Call iCopy(nSym,NumCho,1,NumCho_Bak,1)
      Else
         ip_InfVec_Bak = 0
         l_InfVec_Bak = 0
         Call iZero(NumCho_Bak,nSym)
      End If

      End
