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
* Copyright (C) 2015, Lasse Kragh Soerensen                            *
*               2015, Roland Lindh                                     *
************************************************************************
      SubRoutine OMQInt(Alpha,nAlpha,Beta, nBeta,Zeta,ZInv,rKappa,P,
     &                  Final,nZeta,nIC,nComp,la,lb,A,RB,nHer,
     &                  Array,nArr,Ccoor,nOrdOp,lOper,iChO,
     &                  iStabM,nStabM,
     &                  PtChrg,nGrid,iAddPot)
************************************************************************
*                                                                      *
* Object: kernel routine for the computation of orbital magnetic       *
*         quadrupole integrals => OMQInt                               *
*                                                                      *
* Called from: OneEl                                                   *
*                                                                      *
* Calling    : QEnter                                                  *
*              RecPrt                                                  *
*              MltPrm                                                  *
*              Util2                                                   *
*              DCopy   (ESSL)                                          *
*              GetMem                                                  *
*              QExit                                                   *
*                                                                      *
*     Author: Lasse Kragh Soerensen and Roland Lindh  2015             *
*             Based on OAMInt                                          *
************************************************************************
      Implicit Real*8 (A-H,O-Z)
#include "real.fh"
#include "itmax.fh"
#include "info.fh"
#include "WrkSpc.fh"
#include "print.fh"
      Real*8 Final(nZeta,(la+1)*(la+2)/2,(lb+1)*(lb+2)/2,nIC),
     &       Zeta(nZeta), ZInv(nZeta), Alpha(nAlpha), Beta(nBeta),
     &       rKappa(nZeta), P(nZeta,3), A(3), RB(3), TC(3),
     &       Array(nZeta*nArr), Ccoor(3)
      Integer iStabM(0:nStabM-1), iStabO(0:7), iDCRT(0:7),
     &          lOper(nComp), iChO(nComp)
*
*     Statement function for Cartesian index
*
      nElem(ixyz) = (ixyz+1)*(ixyz+2)/2
*
      iRout = 210
      iPrint = nPrint(iRout)
      Call qEnter('OMQInt')
*
      nip = 1
      ipB = nip
      nip = nip + nZeta

! L + 1 component
      ipS1 = nip
      nip = nip + nZeta*nElem(la)*nElem(lb+1)*6 ! not ncomp

! L - 1 component
      ipS2 =  1
      If (lb.gt.0) Then
         ipS2 = nip
         nip = nip + nZeta*nElem(la)*nElem(lb-1)*6 ! not ncomp
      End If

! L + 0 component
      ipS3 = nip
      nip = nip + nZeta*nElem(la)*nElem(lb)*3

      ipRes = nip
      nip = nip + nZeta*nElem(la)*nElem(lb)*nComp
      If (nip-1.gt.nZeta*nArr) Then
         Call WarningMessage(2,' OMQInt: nip-1.gt.nZeta*nArr')
         Call Abend()
      End If
      ipArr = nip
      mArr = (nArr*nZeta - (nip-1))/nZeta
*
      Call DCopy_(nZeta*nElem(la)*nElem(lb)*nIC,Zero,0,Final,1)
*
      llOper = lOper(1)
      Do 90 iComp = 2, nComp
         llOper = iOr(llOper,lOper(iComp))
 90   Continue
      Call SOS(iStabO,nStabO,llOper)
      Call DCR(LmbdT,iOper,nIrrep,iStabM,nStabM,iStabO,nStabO,
     &         iDCRT,nDCRT)
*
      ipOff = ipB
      Do 100 iAlpha = 1, nAlpha
         Call DCopy_(nBeta,Beta,1,Array(ipOff),nAlpha)
         ipOff = ipOff + 1
 100  Continue
*
      Do 102 lDCRT = 0, nDCRT-1
         TC(1) = DBLE(iPhase(1,iDCRT(lDCRT)))*Ccoor(1)
         TC(2) = DBLE(iPhase(2,iDCRT(lDCRT)))*Ccoor(2)
         TC(3) = DBLE(iPhase(3,iDCRT(lDCRT)))*Ccoor(3)
*
         iComp=6 ! Why are these here ncomp is passed down?
*
         nHer = (la + (lb+1) + (nOrdOp-1) + 2) / 2
         Call MltPrm(Alpha,nAlpha,Beta,nBeta,Zeta,ZInv,rKappa,P,
     &             Array(ipS1),nZeta,iComp,la,lb+1,A,RB,nHer,
     &             Array(ipArr),mArr,TC,nOrdOp-1)
*
         If (lb.gt.0) Then
            nHer = (la + (lb-1) + (nOrdOp-1) + 2) / 2
            Call MltPrm(Alpha,nAlpha,Beta,nBeta,Zeta,ZInv,rKappa,P,
     &                Array(ipS2),nZeta,iComp,la,lb-1,A,RB,nHer,
     &                Array(ipArr),mArr,TC,nOrdOp-1)
         End If
*
         iComp=3 ! Why are these here ncomp is passed down?
*
         nHer = (la + lb + (nOrdOp-2) + 2) / 2
!        check to see dipole integral sure looks a lot like dipole integrals
         Call MltPrm(Alpha,nAlpha,Beta,nBeta,Zeta,ZInv,rKappa,P,
     &             Array(ipS3),nZeta,iComp,la,lb,A,RB,nHer,
     &             Array(ipArr),mArr,TC,nOrdOp-2)
*
*        Combine derivatives of dipole integrals to generate the
*        orbital angular momentum integrals.
*
         Call Util3(Array(ipB),nZeta,Array(ipRes),la,lb,Array(ipS1),
     &                                      Array(ipS3),Array(ipS2))
*
*--------Accumulate contributions
*
         nOp = NrOpr(iDCRT(lDCRT),iOper,nIrrep)
         Call SymAdO(Array(ipRes),nZeta,la,lb,nComp,Final,nIC,
     &               nOp         ,lOper,iChO,One)
*
 102  Continue
*
      Call qExit('OMQInt')
      Return
c Avoid unused argument warnings
      If (.False.) Then
         Call Unused_real(PtChrg)
         Call Unused_integer(nGrid)
         Call Unused_integer(iAddPot)
      End If
      End
