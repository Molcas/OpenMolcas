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
      SubRoutine Cho_XCV_WrVec(irc,Vec,l_Vec,NVT,l_NVT,
     &                         myRankSP,l_myRankSP,SP)
C
C     Thomas Bondo Pedersen, April 2010.
C
C     Purpose: write partial Cholesky vectors to disk.
C              (Parallel two-step algorithm)
C
#if defined (_DEBUGPRINT_)
      use ChoSwp, only: nnBstRSh
#endif
      Implicit None
      Integer irc
      Integer l_Vec, l_NVT, l_myRankSP
      Real*8  Vec(l_Vec)
      Integer NVT(l_NVT)
      Integer myRankSP(l_myRankSP)
      Integer SP
#include "cho_para_info.fh"

#if defined (_DEBUGPRINT_)
#include "cholesky.fh"
#include "choptr.fh"
#include "WrkSpc.fh"

      Integer iSym, iSP, n

      If (l_NVT.lt.nSym .or. l_myRankSP.lt.nnShl) Then
         irc=-1
         Return
      End If
      If (SP.lt.1 .or. SP.gt.nnShl) Then
         irc=-2
         Return
      End If
      iSP=myRankSP(SP)
      n=nnBstRSh(1,iSP,2)*NVT(1)
      Do iSym=2,nSym
         n=n+nnBstRSh(iSym,iSP,2)*NVT(iSym)
      End Do
      If (l_Vec.lt.n) Then
         irc=-3
         Return
      End If

      If (.not.Cho_Real_Par) Then
         n=0
         Do iSym=1,nSym
            If (NVT(iSym).ne.NumCho(iSym)) Then
               n=n+1
            End If
         End Do
         If (n.ne.0) Then
            irc=-4
            Return
         End If
         If (SP.ne.myRankSP(SP)) Then
            irc=-5
            Return
         End If
      End If
#endif

      If (Cho_Real_Par) Then
         ! Parallel: block write to temp files
         Call Cho_XCV_WrVec_Par(irc,Vec,NVT,myRankSP,SP)
      Else
         ! Serial: write directly to vector files
         Call Cho_XCV_WrVec_Ser(irc,Vec,SP)
      End If

      End
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C SERIAL VERSION
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SubRoutine Cho_XCV_WrVec_Ser(irc,Vec,iSP)
C
C     Simply write the partial vectors to disk at the appropriate
C     addresses on the vector files.
C
      use ChoSwp, only: nnBstRSh, iiBstRSh
      Implicit None
      Integer irc
      Real*8  Vec(*)
      Integer iSP
#include "cholesky.fh"
#include "choptr.fh"
#include "WrkSpc.fh"

      Integer iOpt
      Parameter (iOpt=1)

      Integer iSym, kV
      Integer lTot, iAdr, iAdr0, j

      irc=0

      kV=1
      Do iSym=1,nSym
         lTot=nnBstRSh(iSym,iSP,2)
         If (lTot .gt. 0) Then
            iAdr0=iiBstRSh(iSym,iSP,2)
            Do J=1,NumCho(iSym)
               iAdr=iAdr0+nnBstR(iSym,2)*(J-1)
               Call DDAFile(LuCho(iSym),iOpt,Vec(kV),lTot,iAdr)
               kV=kV+lTot
            End Do
         End If
      End Do

      End
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C PARALLEL VERSION
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SubRoutine Cho_XCV_WrVec_Par(irc,Vec,NVT,myRankSP,SP)
C
C     Write the vectors in blocks.
C
      use ChoSwp, only: nnBstRSh
      Implicit None
      Integer irc
      Real*8  Vec(*)
      Integer NVT(*)
      Integer myRankSP(*)
      Integer SP
#include "cholesky.fh"
#include "choptr.fh"
#include "WrkSpc.fh"

      Integer, Parameter:: iOpt=1

      Integer iSym, kV, j
      Integer lTot, iAdr
      Integer iSP

      irc=0

      iSP=myRankSP(SP)
      kV=1
      Do iSym=1,nSym
         lTot=nnBstRSh(iSym,iSP,2)*NVT(iSym)
         If (lTot .gt. 0) Then
            iAdr=0
            Do j=1,SP-1
               iAdr=iAdr+nnBstRSh(iSym,myRankSP(j),2)*NVT(iSym)
            End Do
            Call DDAFile(LuTmp(iSym),iOpt,Vec(kV),lTot,iAdr)
            kV=kV+lTot
         End If
      End Do

      End
