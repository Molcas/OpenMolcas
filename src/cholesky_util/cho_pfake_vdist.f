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
      Use Para_Info, Only: nProcs, Is_Real_Par
      use ChoSwp, only: InfVec
      Implicit None
#include "cholesky.fh"
#include "choglob.fh"
#include "stdalloc.fh"

      Character(LEN=15), Parameter:: SecNam = 'Cho_PFake_VDist'

      Integer iSym
      Integer l_Wrk
      Integer iRedC, iV, nV
      Integer nRead

      Integer, Allocatable:: InfV(:,:), IDV(:)
      Real*8, Allocatable:: Wrk(:)

C     Return if nothing to do.
C     ------------------------

      If (nProcs.eq.1 .or. .not.Is_Real_Par()) Return ! serial run
      If (.not.CHO_FAKE_PAR) Return ! true parallel run

C     Memory allocation.
C     ------------------

      Call mma_allocate(InfV,2,MaxVec+1,Label='InfV')
      Call mma_allocate(IDV,MaxVec,Label='IDV')
      Call mma_maxDBLE(l_Wrk)
      Call mma_allocate(Wrk,l_Wrk,Label='Wrk')

C     Distribute vectors in each symmetry:
C     read from LuCho files; put into LuCho_G files.
C     ----------------------------------------------

      iRedC = -1
      Do iSym = 1,nSym
         InfV(:,:)=0
         nV = 0
         Call Cho_Distrib_Vec(1,NumCho(iSym),IDV,nV)
         iV = 0
         Do While (iV .lt. nV)
            nRead = 0
            Call Cho_PFake_GetVec(Wrk,SIZE(Wrk),
     &                            IDV(iV+1),nV-iV,
     &                            InfV(:,iV+1),
     &                            iSym,nRead,iRedC)
            If (nRead .lt. 1) Then
               Call Cho_Quit('Insufficient memory in '//SecNam,101)
            End If
            Call Cho_PFake_PutVec(Wrk,InfV,nRead,iSym,iV+1)
            iV = iV + nRead
         End Do
#if defined (_DEBUGPRINT_)
         If (iV .ne. nV) Then
            Call Cho_Quit('Logical error in '//SecNam,103)
         End If
#endif
         Do iV = 1,nV
            InfVec(iV,3,iSym)= InfV(2,iV)
         End Do
      End Do

C     Deallocation.
C     -------------

      Call mma_deallocate(Wrk)
      Call mma_deallocate(IDV)
      Call mma_deallocate(InfV)

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
