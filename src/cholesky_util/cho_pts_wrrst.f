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
      SubRoutine Cho_PTS_WrRst(irc,NVT,l_NVT)
C
C     Thomas Bondo Pedersen, April 2010.
C
C     Purpose: Write restart files (parallel two-step algorithm).
C
      use ChoSwp, only: InfRed, InfVec
      Implicit None
      Integer irc
      Integer l_NVT
      Integer NVT(l_NVT)
#include "cholesky.fh"
#include "WrkSpc.fh"

      Integer N2

#if defined (_DEBUGPRINT_)
      Character*13 SecNam
      Parameter (SecNam='Cho_PTS_WrRst')
      Integer ip_IDV, l_IDV
      Integer myNumCho(8)
#endif

      Integer iSym
      Integer, Pointer:: InfVcT(:,:,:)
      Integer iV, iAdr

*                                                                      *
************************************************************************
*                                                                      *
      Interface
      Subroutine Cho_X_GetIP_InfVec(InfVcT)
      Integer, Pointer:: InfVct(:,:,:)
      End Subroutine Cho_X_GetIP_InfVec
      End Interface
*                                                                      *
************************************************************************
*                                                                      *
      ! Init return code
      irc=0
      N2 = SIZE(InfVec,2)

#if defined (_DEBUGPRINT_)
      ! check that NumCho agrees with distribution
      If (l_NVT.lt.nSym) Then
         irc=-1
         Return
      End If
      Do iSym=1,nSym
         l_IDV=NVT(iSym)
         Call GetMem('IDV','Allo','Inte',ip_IDV,l_IDV)
         myNumCho(iSym)=0
         Call Cho_P_Distrib_Vec(1,NVT(iSym),iWork(ip_IDV),
     &                          myNumCho(iSym))
         Call GetMem('IDV','Free','Inte',ip_IDV,l_IDV)
         If (NumCho(iSym) .ne. myNumCho(iSym)) Then
            Write(LuPri,*)
     &      SecNam,': NumCho discrepancy in sym. ',iSym
            Write(LuPri,*)
     &      '  NumCho=',NumCho(iSym)
            Write(LuPri,*)
     &      'myNumCho=',myNumCho(iSym)
            Write(LuPri,*)
     &      '     NVT=',NVT(iSym)
            irc=1
         End If
      End Do
      Call iCopy(nSym,NumCho,1,myNumCho,1)
      Call Cho_GAIGOp(myNumCho,nSym,'+')
      Do iSym=1,nSym
         If (myNumCho(iSym).ne.NVT(iSym)) Then
            Write(LuPri,*)
     &      SecNam,': NumCho discrepancy in sym. ',iSym
            Write(LuPri,*)
     &      'Sum of NumCho across nodes=',myNumCho(iSym)
            Write(LuPri,*)
     &      '                       NVT=',NVT(iSym)
            irc=2
         End If
      End Do
      If (irc.ne.0) Return
#endif

      ! Erase existing files from first step
      ! (Keep map and vector files, though)
      If (LuRst.lt.1) Then
         LuRst=7
         Call DAName_MF_WA(LuRst,'CHORST')
      End If
      Call DAEras(LuRst)
      LuRst=0
      If (LuRed.lt.1) Then
         LuRed=7
         Call DAName_MF_WA(LuRst,'CHRED')
      End If
      Call DAEras(LuRed)
      LuRed=0
      If (LuMap.gt.0) Then
         Call DAClos(LuMap)
         LuMap=0
      End If
      Do iSym=1,nSym
         If (LuCho(iSym).gt.0) Then
            Call DAClos(LuCho(iSym))
            LuCho(iSym)=0
         End If
      End Do

      ! Re-open files
      Call Cho_OpenVR(1,2)

      ! Set InfRed corresponding to only one reduced set
      ! (All vectors are now stored in 1st reduced set)
      InfRed(1)=0
      ! Set InfVec data (only disk addresses need updating)
      Call Cho_X_GetIP_InfVec(InfVcT)
      Do iSym=1,nSym
         ! InfVec(iV,1,iSym): parent diagonal in rs1 of vector iV of
         ! symmetry iSym
         Do iV=1,NVT(iSym)
            InfVec(iV,1,iSym) = InfVcT(iV,1,iSym)
         End Do
         Do iV=NVT(iSym)+1,MaxVec
            InfVec(iV,1,iSym)=0
         End Do
         ! InfVec(iV,2,iSym): reduced set of vector iV of
         ! symmetry iSym (always rs1 in this implementation)
         Do iV=1,NVT(iSym)
            InfVec(iV,2,iSym)=1
         End Do
         Do iV=NVT(iSym)+1,MaxVec
            InfVec(iV,2,iSym)=0
         End Do
         ! InfVec(iV,3,iSym): disk address of vector iV of
         ! symmetry iSym (always rs1 in this implementation)
         ! Note: this information is for the vectors on this
         ! node ONLY!!
         iAdr=0
         Do iV=1,NumCho(iSym)
            InfVec(iV,3,iSym)=iAdr
            iAdr=iAdr+nnBstR(iSym,1)
         End Do
         Do iV=NumCho(iSym)+1,MaxVec
            InfVec(iV,3,iSym)=-1
         End Do
         ! InfVec(iV,4,iSym): WA disk address of vector iV of
         ! symmetry iSym (redundant? Oh yes, but used for statistics...)
         iAdr=0
         Do iV=1,NVT(iSym)
            InfVec(iV,4,iSym)=iAdr
            iAdr=iAdr+nnBstR(iSym,1)
         End Do
         Do iV=NVT(iSym)+1,MaxVec
            InfVec(iV,4,iSym)=-1
         End Do
         ! InfVec(iV,5,iSym): is defined automatically by Cho_X_Init
         ! No need to set it here - use a dummy value
         InfVec(:,5,iSym)=-1
      End Do

      ! Make sure that all vector info is written
      Call iSwap(nSym,NVT,1,NumCho,1)

      ! Write restart file
      Call Cho_WrRstC(1)

      ! Write reduced set
      Call Cho_PutRed(1,1)

      ! Swap back
      Call iSwap(nSym,NVT,1,NumCho,1)

      End
