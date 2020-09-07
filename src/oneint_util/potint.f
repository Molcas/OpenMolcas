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
      SubRoutine PotInt(Alpha,nAlpha,Beta, nBeta,Zeta,ZInv,rKappa,P,
     &                  Final,nZeta,nIC,nComp,la,lb,A,RB,nRys,
     &                  Array,nArr,CCoor,nOrdOp,lOper,iChO,
     &                  iStabM,nStabM,PtChrg,nGrid,iAddPot)
************************************************************************
*                                                                      *
* Object: kernel routine for the computation of potential integrals    *
*                                                                      *
* Called from: OneEl                                                   *
*                                                                      *
* Calling    : QEnter                                                  *
*              RecPrt                                                  *
*              DCopy   (ESSL)                                          *
*              mHrr                                                    *
*              DCR                                                     *
*              Rys                                                     *
*              Hrr                                                     *
*              DaXpY   (ESSL)                                          *
*              QExit                                                   *
*                                                                      *
*     Author: Roland Lindh, Dept. of Theoretical Chemistry, University *
*             of Lund, Sweden, January '91                             *
************************************************************************
      Implicit Real*8 (A-H,O-Z)
*     Used for normal nuclear attraction integrals
      External TNAI, Fake, XCff2D, XRys2D
#include "itmax.fh"
#include "info.fh"
#include "real.fh"
#include "oneswi.fh"
#include "print.fh"
      Real*8 Final(nZeta,(la+1)*(la+2)/2,(lb+1)*(lb+2)/2,nIC),
     &       Zeta(nZeta), ZInv(nZeta), Alpha(nAlpha), Beta(nBeta),
     &       rKappa(nZeta), P(nZeta,3), A(3), RB(3), CCoor(3,*),
     &       Array(nZeta*nArr),ptchrg(nGrid)
      Integer lOper(nComp), iStabM(0:nStabM-1), iStabO(0:7), iPh(3)
*-----Local arrys
      Real*8 TC(3), Coora(3,4), Coori(3,4), CoorAC(3,2)
      Logical EQ, NoSpecial
      Integer iAnga(4), iDCRT(0:7), iChO(nComp)
*
*     Statement function for Cartesian index
*
      nElem(ixyz) = (ixyz+1)*(ixyz+2)/2
      nabSz(ixyz) = (ixyz+1)*(ixyz+2)*(ixyz+3)/6  - 1
*
      Call fzero(final,nZeta*nElem(la)*nElem(lb)*nIC)
      len=nZeta*nElem(la)*nElem(lb)*nIC
*
      iAnga(1) = la
      iAnga(2) = lb
      iAnga(3) = 0
      iAnga(4) = 0
      call dcopy_(3,A,1,Coora(1,1),1)
      call dcopy_(3,RB,1,Coora(1,2),1)
      call dcopy_(2*3,Coora,1,Coori,1)
      mabMin = nabSz(Max(la,lb)-1)+1
      mabMax = nabSz(la+lb)
      If (EQ(A,RB)) mabMin=nabSz(la+lb-1)+1
*
*     Compute FLOP's and size of work array which Hrr will use.
*
      Call mHrr(la,lb,nFLOP,nMem)
*
*     Find center to accumulate angular momentum on. (HRR)
*
      If (la.ge.lb) Then
         call dcopy_(3,A,1,CoorAC(1,1),1)
      Else
         call dcopy_(3,RB,1,CoorAC(1,1),1)
      End If
*
      llOper = lOper(1)
*
*     Loop over grid
*
*-----------Find the DCR for M and S
*
      Call SOS(iStabO,nStabO,llOper)
      Call DCR(LmbdT,iOper,nIrrep,iStabM,nStabM,iStabO,nStabO,
     &         iDCRT,nDCRT)
c     Fact = DBLE(nStabM) / DBLE(LmbdT)
      FACT=1.D0

      nT = nZeta
      Chrg=-1.d0
      NoSpecial=.True.
*
      Do lDCRT = 0, nDCRT-1

         Do i = 1, 3
            iph(i) = iPhase(i,iDCRT(lDCRT))
         End Do
         nOp = NrOpr(iDCRT(lDCRT),iOper,nIrrep)

         Do 100 iGrid = 1, nGrid
            If (iAddPot.ne.0) Chrg=ptchrg(iGrid)
            If (Chrg.eq.Zero) Go To 100
*
               Do i = 1, 3
                  TC(i)=DBLE(iPh(i))*CCoor(i,iGrid)
                  CoorAC(i,2)=TC(i)
                  Coori(i,3)=TC(i)
                  Coori(i,4)=TC(i)
                  Coora(i,3)=TC(i)
                  Coora(i,4)=TC(i)
               End Do
*
*              Compute integrals with the Rys quadrature.
*
               Call Rys(iAnga,nT,Zeta,ZInv,nZeta,
     &                  [One],[One],1,P,nZeta,
     &                  TC,1,rKappa,[One],Coori,Coora,CoorAC,
     &                  mabmin,mabmax,0,0,Array,nArr*nZeta,
     &                  TNAI,Fake,XCff2D,XRys2D,NoSpecial)
*
*--------------Use the HRR to compute the required primitive integrals.
*
               Call HRR(la,lb,A,RB,Array,nZeta,nMem,ipIn)
*
*--------------Accumulate contributions to the symmetry adapted operator
*
               Call SymAdO(Array(ipIn),nZeta,la,lb,nComp,Final,nIC,
     &                     nOp,lOper,iChO,-Fact*Chrg)
 100        Continue
      End Do
*
      Return
c Avoid unused argument warnings
      If (.False.) Then
         Call Unused_real_array(Alpha)
         Call Unused_real_array(Beta)
         Call Unused_integer(nRys)
         Call Unused_integer(nOrdOp)
      End If
      End
      SubRoutine Pot_nuc(CCoor,pot,nGrid)
      use Basis_Info
      Implicit Real*8 (A-H,O-Z)
#include "itmax.fh"
#include "info.fh"
#include "real.fh"
#include "print.fh"
      Real*8  CCoor(3,nGrid),pot(nGrid)
      Real*8 C(3), TC(3)
      Integer iStabM(0:7),iDCRT(0:7)
*
*  compute nuclear contribution to potential
*
      kdc = 0
      Do iGrid=1,nGrid
         pot(iGrid)=0d0
      End Do
*
chjw is this always correct?
      istabm(0)=0
      nstabm=1
*
      Do 100 kCnttp = 1, nCnttp
         If (dbsc(kCnttp)%Charge.eq.Zero) Go To 111
*
         Do 101 kCnt = 1, dbsc(kCnttp)%nCntr
*
            C(1:3) = dbsc(kCnttp)%Coor(1:3,kCnt)
            Call DCR(LmbdT,iOper,nIrrep,iStabM,nStabM,
     &               jStab(0,kdc+kCnt) ,nStab(kdc+kCnt),iDCRT,nDCRT)
            Fact = DBLE(nStabM) / DBLE(LmbdT)
*
            Do lDCRT = 0, nDCRT-1
               TC(1) = DBLE(iPhase(1,iDCRT(lDCRT)))*C(1)
               TC(2) = DBLE(iPhase(2,iDCRT(lDCRT)))*C(2)
               TC(3) = DBLE(iPhase(3,iDCRT(lDCRT)))*C(3)
*
               Do iGrid=1,nGrid
                 r12=sqrt((TC(1)-CCoor(1,iGrid))**2
     &                   +(TC(2)-CCoor(2,iGrid))**2
     &                   +(TC(3)-CCoor(3,iGrid))**2)
                 if(r12.gt.1.d-8)
     &            pot(iGrid)=pot(iGrid)+dbsc(kCnttp)%Charge*fact/r12
               End Do
*
            End Do
 101     Continue
 111     kdc = kdc + dbsc(kCnttp)%nCntr
 100  Continue
*
      Return
      End
