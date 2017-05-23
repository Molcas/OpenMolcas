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
* Copyright (C) 2007, Thomas Bondo Pedersen                            *
************************************************************************
      SubRoutine Cho_PFake_VDist()
C
C     Author: Thomas Bondo Pedersen, April 2007.
C     Purpose: fake parallel distribution of vectors.
C
      Implicit None
#include "cholesky.fh"
#include "choptr.fh"
#include "choglob.fh"
#include "WrkSpc.fh"
#include "cho_para_info.fh"

      Character*15 SecNam
      Parameter (SecNam = 'Cho_PFake_VDist')

      Integer iSym
      Integer ip_IDV, l_IDV
      Integer ip_Wrk, l_Wrk
      Integer ip_InfV, l_InfV
      Integer iRedC, iV, nV
      Integer kOff, nRead

      Integer i, j
      Integer InfV
      InfV(i,j)=iWork(ip_InfV-1+2*(j-1)+i)

C     Return if nothing to do.
C     ------------------------

      If (nProcs.eq.1 .or. .not.Is_Real_Par()) Return ! serial run
      If (.not.CHO_FAKE_PAR) Return ! true parallel run

C     Memory allocation.
C     ------------------

      l_InfV = 2*(MaxVec+1)
      l_IDV = MaxVec
      Call GetMem('InfV_','Allo','Inte',ip_InfV,l_InfV)
      Call GetMem('IDV_','Allo','Inte',ip_IDV,l_IDV)
      Call GetMem('MAX_','Max ','Real',ip_Wrk,l_Wrk)
      Call GetMem('Wrk_','Allo','Real',ip_Wrk,l_Wrk)

C     Distribute vectors in each symmetry:
C     read from LuCho files; put into LuCho_G files.
C     ----------------------------------------------

      iRedC = -1
      Do iSym = 1,nSym
         Call iZero(iWork(ip_InfV),l_InfV)
         nV = 0
         Call Cho_Distrib_Vec(1,NumCho(iSym),iWork(ip_IDV),nV)
         iV = 0
         Do While (iV .lt. nV)
            nRead = 0
            Call Cho_PFake_GetVec(Work(ip_Wrk),l_Wrk,
     &                            iWork(ip_IDV+iV),nV-iV,
     &                            iWork(ip_InfV+2*iV),
     &                            iSym,nRead,iRedC)
            If (nRead .lt. 1) Then
               Call Cho_Quit('Insufficient memory in '//SecNam,101)
            End If
            Call Cho_PFake_PutVec(Work(ip_Wrk),iWork(ip_InfV),nRead,
     &                            iSym,iV+1)
            iV = iV + nRead
         End Do
#if defined (_DEBUG_)
         If (iV .ne. nV) Then
            Call Cho_Quit('Logical error in '//SecNam,103)
         End If
#endif
         kOff = ip_InfVec - 1 + MaxVec*InfVec_N2*(iSym-1) + MaxVec*2
         Do iV = 1,nV
            iWork(kOff+iV) = InfV(2,iV)
         End Do
      End Do

C     Deallocation.
C     -------------

      Call GetMem('Wrk_','Free','Real',ip_Wrk,l_Wrk)
      Call GetMem('IDV_','Free','Inte',ip_IDV,l_IDV)
      Call GetMem('InfV_','Free','Inte',ip_InfV,l_InfV)

      End
      SubRoutine Cho_PFake_GetVec(Vec,lVec,IDV,lIDV,InfV,iSym,nRead,
     &                            iRedC)
      Implicit None
      Integer lVec, lIDV
      Real*8  Vec(lVec)
      Integer IDV(lIDV)
      Integer InfV(2,*)
      Integer iSym, nRead, iRedC

      Character*16 SecNam
      Parameter (SecNam = 'Cho_PFake_GetVec')

      Integer ipV, Mem, iVec, n, m

      nRead = 0
      ipV = 1
      Mem = lVec
      Do iVec = 1,lIDV
         n = 0
         m = 0
         Call Cho_VecRd(Vec(ipV),Mem,IDV(iVec),IDV(iVec),iSym,
     &                  n,iRedC,m)
         If (n .eq. 1) Then
            nRead = nRead + 1
            ipV = ipV + m
            Mem = Mem - m
            InfV(1,iVec) = m
         Else If (n .eq. 0) Then
            Return
         Else
            Call Cho_Quit('Logical error in '//SecNam,103)
         End If
      End Do

      End
      SubRoutine Cho_PFake_PutVec(Vec,InfV,nVec,iSym,iV1)
      Implicit None
      Real*8 Vec(*)
      Integer InfV(2,*)
      Integer nVec, iSym, iV1
#include "cholesky.fh"
#include "choglob.fh"

      Character*16 SecNam
      Parameter (SecNam = 'Cho_PFake_PutVec')

      Integer iOpt, lTot, iAdr
      Integer iVec, iPos

      If (nVec .lt. 1) Return

      If (CHO_ADRVEC .eq. 1) Then
         iOpt = 1
         lTot = InfV(1,iV1)
         Do iVec = iV1+1,iV1+nVec-1
            lTot = lTot + InfV(1,iVec)
         End Do
         iAdr = InfV(2,iV1)
         Call dDAFile(LuCho_G(iSym),iOpt,Vec,lTot,iAdr)
         Do iVec = iV1+1,iV1+nVec
            InfV(2,iVec) = InfV(2,iVec-1) + InfV(1,iVec-1)
         End Do
      Else If (CHO_ADRVEC .eq. 2) Then
         iPos = 1
         Do iVec = iV1,iV1+nVec-1
            iOpt = 1
            lTot = InfV(1,iVec)
            iAdr = InfV(2,iVec)
            Call dDAFile(LuCho_G(iSym),iOpt,Vec(iPos),lTot,iAdr)
            iPos = iPos + InfV(1,iVec)
            InfV(2,iVec+1) = iAdr
         End Do
      Else
         Call Cho_Quit('Illegal CHO_ADRVEC in '//SecNam,102)
      End If

      End
