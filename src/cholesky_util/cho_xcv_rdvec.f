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
      SubRoutine Cho_XCV_RdVec(irc,Vec,l_Vec,NVT,myRankSP,n_myRankSP,
     &                         J1,J2,iSym)
C
C     Thomas Bondo Pedersen, April 2010.
C
C     Purpose: Read partial Cholesky vectors J1 to J2 on disk.
C              (Parallel two-step algorithm)
C
#if defined (_DEBUGPRINT_)
      use ChoSwp, only: nnBstRSh
#endif
      Implicit None
      Integer irc
      Integer l_Vec
      Real*8  Vec(l_Vec)
      Integer NVT
      Integer n_myRankSP
      Integer myRankSP(n_myRankSP)
      Integer J1, J2, iSym

#if defined (_DEBUGPRINT_)
#include "cho_para_info.fh"
#include "cholesky.fh"
#include "choptr.fh"
#include "WrkSpc.fh"

      Integer i, n
#endif

      irc=0
      If (n_myRankSP.eq.0 .or. (J2-J1+1).eq.0) Then
         Return ! nothing to do
      End If

#if defined (_DEBUGPRINT_)
      If (n_myRankSP.lt.1 .or. n_myRankSP.gt.nnShl .or. iSym.lt.1 .or.
     &    iSym.gt.nSym .or. NVT.lt.1) Then
         irc=-1
         Return
      End If
      If (J1.lt.1 .or. J1.gt.NVT .or. J2.lt.1 .or. J2.gt.NVT .or.
     &    J1.gt.J2 .or. (J2-J1+1).gt.NVT) Then
         irc=-2
         Return
      End If
      n=nnBstRSh(iSym,myRankSP(1),2)*(J2-J1+1)
      Do i=2,n_myRankSP
         n=n+nnBstRSh(iSym,myRankSP(i),2)*(J2-J1+1)
      End Do
      If (l_Vec.lt.n) Then
         irc=-3
         Return
      End If

      If (.not.Cho_Real_Par) Then
         If (NVT.ne.NumCho(iSym)) Then
            irc=-4
            Return
         End If
         n=0
         Do i=1,n_myRankSP
            If (myRankSP(i).ne.i) Then
               n=n+1
            End If
         End Do
         If (n.ne.0) Then
            irc=-5
            Return
         End If
      End If
#endif

      ! Block read on temp files
      Call Cho_XCV_RdVec_(irc,Vec,myRankSP,n_myRankSP,NVT,J1,J2,iSym)

      End
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C BLOCK READ
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SubRoutine Cho_XCV_RdVec_(irc,Vec,myRankSP,n_myRankSP,NVT,
     &                          J1,J2,iSym)
C
C     Read the vector blocks.
C
      use ChoSwp, only: nnBstRSh
      Implicit None
      Integer irc
      Real*8  Vec(*)
      Integer n_myRankSP
      Integer myRankSP(n_myRankSP)
      Integer NVT
      Integer J1, J2, iSym
#include "cholesky.fh"
#include "choptr.fh"
#include "WrkSpc.fh"

      Integer iOpt
      Parameter (iOpt=2)

      Integer kV, n, i
      Integer lTot, iAdr, iAdr0
      Integer iSP

      irc=0

      n=J2-J1+1
      iAdr0=0
      kV=1
      Do i=1,n_myRankSP
         iSP=myRankSP(i)
         lTot=nnBstRSh(iSym,iSP,2)*n
         If (lTot.gt.0) Then
            iAdr=iAdr0+nnBstRSh(iSym,iSP,2)*(J1-1)
            Call DDAFile(LuTmp(iSym),iOpt,Vec(kV),lTot,iAdr)
            kV=kV+lTot
         End If
         iAdr0=iAdr0+nnBstRSh(iSym,iSP,2)*NVT
      End Do

      End
