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
* Copyright (C) 1991, Roland Lindh                                     *
************************************************************************
      SubRoutine EPEInt(Alpha,nAlpha,Beta, nBeta,Zeta,ZInv,rKappa,P,
     &                  Final,nZeta,nIC,nComp,la,lb,A,RB,nRys,
     &                  Array,nArr,Ccoor,nOrdOp,lOper,iChO,
     &                  iStabM,nStabM,
     &                  PtChrg,nGrid,iAddPot)
************************************************************************
*                                                                      *
* Object: kernel routine for the computation of nuclear attraction     *
*         integrals.                                                   *
*                                                                      *
* Called from: EFInt                                                   *
*                                                                      *
* Calling    : QEnter                                                  *
*              RecPrt                                                  *
*              DCopy  (ESSL)                                           *
*              ICopy                                                   *
*              SOS                                                     *
*              Rys                                                     *
*              Hrr                                                     *
*              QExit                                                   *
*                                                                      *
*     Author: Roland Lindh, Dept. of Theoretical Chemistry, University *
*             of Lund, Sweden, February '91                            *
************************************************************************
      Implicit Real*8 (A-H,O-Z)
      External TNAI, Fake, Cff2D, XRys2D
#include "real.fh"
#include "itmax.fh"
#include "info.fh"
#include "WrkSpc.fh"
#include "print.fh"
      Real*8 Final(nZeta,(la+1)*(la+2)/2,(lb+1)*(lb+2)/2,nIC),
     &       Zeta(nZeta), ZInv(nZeta), Alpha(nAlpha), Beta(nBeta),
     &       rKappa(nZeta), P(nZeta,3), A(3), RB(3), TC(3),
     &       Array(nZeta*nArr), Ccoor(3)
      Real*8 Coori(3,4), Coora(3,4), CoorAC(3,2)
      Integer iAnga(4), iStabM(0:nStabM-1), iDCRT(0:7),
     &          iStabO(0:7), iChO(nComp), lOper(nComp)
      Logical EQ, NoSpecial
*
*     Statement function for Cartesian index
*
      nElem(ixyz) = (ixyz+1)*(ixyz+2)/2
      nabSz(ixyz) = (ixyz+1)*(ixyz+2)*(ixyz+3)/6  - 1
*
      iRout = 201
      iPrint = nPrint(iRout)
      Call qEnter('EPEInt')
*
      call dcopy_(nZeta*nElem(la)*nElem(lb)*nIC,Zero,0,Final,1)
*
      iAnga(1) = la
      iAnga(2) = lb
      iAnga(3) = 0
      iAnga(4) = 0
      call dcopy_(3,A,1,Coora(1,1),1)
      call dcopy_(3,RB,1,Coora(1,2),1)
      call dcopy_(2*3,Coora(1,1),1,Coori(1,1),1)
      mabMin = nabSz(Max(la,lb)-1)+1
      mabMax = nabSz(la+lb)
      If (EQ(A,RB)) mabMin = nabSz(la+lb-1)+1
*
*     Compute FLOP's and size of the work array which Hrr will use.
*
      Call mHrr(la,lb,nFLOP,nMem)
*
*-----Find center to accumulate angular momentum on.
*
      If (la.ge.lb) Then
         call dcopy_(3,A,1,CoorAC(1,1),1)
      Else
         call dcopy_(3,RB,1,CoorAC(1,1),1)
      End If
*
      llOper = lOper(1)
      Do 90 iComp = 2, nComp
         llOper = iOr(llOper,lOper(iComp))
 90   Continue
      Call SOS(iStabO,nStabO,llOper)
      Call DCR(LmbdT,iOper,nIrrep,iStabM,nStabM,iStabO,nStabO,
     &         iDCRT,nDCRT)
      Do 100 lDCRT = 0, nDCRT-1
         TC(1) = DBLE(iPhase(1,iDCRT(lDCRT)))*CCoor(1)
         TC(2) = DBLE(iPhase(2,iDCRT(lDCRT)))*CCoor(2)
         TC(3) = DBLE(iPhase(3,iDCRT(lDCRT)))*CCoor(3)
         call dcopy_(3,TC,1,CoorAC(1,2),1)
         call dcopy_(3,TC,1,Coori(1,3),1)
         call dcopy_(3,TC,1,Coori(1,4),1)
         call dcopy_(3,TC,1,Coora(1,3),1)
         call dcopy_(3,TC,1,Coora(1,4),1)
*
*--------Compute primitive integrals before the application of HRR.
*
         nT = nZeta
         NoSpecial=.True.
         Call Rys(iAnga,nt,Zeta,ZInv,nZeta,One,One,1,P,nZeta,
     &            TC,1,rKappa,One,Coori,Coora,CoorAC,
     &            mabmin,mabmax,0,0,Array,nArr*nZeta,
     &            TNAI,Fake,Cff2D,xRys2D,NoSpecial)
*
*--------Use the HRR to compute the required primitive integrals.
*
         Call HRR(la,lb,A,RB,Array,nZeta,nMem,ipIn)
*
*--------Accumulate contributions to the symmetry adaped operator
*
         nOp = NrOpr(iDCRT(lDCRT),iOper,nIrrep)
         Call SymAdO(Array(ipIn),nZeta,la,lb,nComp,Final,nIC,
     &               nOp         ,lOper,iChO,One)
*
 100  Continue
*
*     Call GetMem(' Exit EPEInt','LIST','REAL',iDum,iDum)
      Call qExit('EPEInt')
      Return
c Avoid unused argument warnings
      If (.False.) Then
         Call Unused_real_array(Alpha)
         Call Unused_real_array(Beta)
         Call Unused_integer(nRys)
         Call Unused_integer(nOrdOp)
         Call Unused_real(PtChrg)
         Call Unused_integer(nGrid)
         Call Unused_integer(iAddPot)
      End If
      End
