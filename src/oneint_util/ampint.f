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
* Copyright (C) 1996, Per Ake Malmqvist                                *
*               1996, Roland Lindh                                     *
************************************************************************
      SubRoutine AMPInt(Alpha,nAlpha,Beta, nBeta,Zeta,ZInv,rKappa,P,
     &                  Final,nZeta,nIC,nComp,la,lb,A,RB,nHer,
     &                  Array,nArr,Ccoor,nOrdOp,lOper,iChO,
     &                  iStabM,nStabM,
     &                  PtChrg,nGrid,iAddPot)
************************************************************************
*                                                                      *
* Object: kernel routine for computing matrix elements of the          *
*         six hermitized products of two angular momentum ops          *
*                                                                      *
* Called from: OneEl                                                   *
*                                                                      *
* Calling    : QEnter                                                  *
*              RecPrt                                                  *
*              MltPrm                                                  *
*              AMPr                                                    *
*              DCopy   (ESSL)                                          *
*              GetMem                                                  *
*              QExit                                                   *
*                                                                      *
*     Author: Per-Ake Malmqvist, Dept. of Theoretical Chemistry,       *
*             University of Lund, SWEDEN                               *
*             November '96                                             *
*     After pattern of other SEWARD soubroutines by R. Lindh.          *
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
      nElem(ixyz) = ((ixyz+1)*(ixyz+2))/2
*
      iRout = 220
      iPrint = nPrint(iRout)
      Call qEnter('AMPInt')
*
      nip = 1
      ipB = nip
      nip = nip + nZeta
      ipTpp = nip
      nip = nip + nZeta*nElem(la)*nElem(lb+2)*6
      ipTp  = nip
      nip = nip + nZeta*nElem(la)*nElem(lb+1)*3
      ipT   = nip
      nip = nip + nZeta*nElem(la)*nElem(lb  )*6
      ipTm  = 1
      ipTmm = 1
      if(lb.gt.0) then
        ipTm  = nip
        nip = nip + nZeta*nElem(la)*nElem(lb-1)*3
        if(lb.gt.1) then
          ipTmm = nip
          nip = nip + nZeta*nElem(la)*nElem(lb-2)*6
        end if
      end if
      ipRes=nip
      nip = nip + nZeta*nElem(la)*nElem(lb)*nComp
      If (nip-1.gt.nZeta*nArr) Then
         Call WarningMessage(2,' AMPInt: nip-1.gt.nZeta*nArr')
         call Abend()
      End If
      ipArr = nip
      mArr = (nArr*nZeta - (nip-1))/nZeta

      call dcopy_(nZeta*nElem(la)*nElem(lb)*nIC,Zero,0,Final,1)

      ipOff = ipB
      Do iAlpha = 1, nAlpha
         call dcopy_(nBeta,Beta,1,Array(ipOff),nAlpha)
         ipOff = ipOff + 1
      End Do

      llOper = lOper(1)
      Do iComp = 2, nComp
         iDum=lOper(iComp)
         llOper = iOr(llOper,iDum)
      End Do

C Compute stabilizer, and then the double coset representation:
      Call SOS(iStabO,nStabO,llOper)
      Call DCR(LmbdT,iOper,nIrrep,iStabM,nStabM,iStabO,nStabO,
     &         iDCRT,nDCRT)

C Loop over the cosets of the stabilizer group:
      Do lDCRT = 0, nDCRT-1
         TC(1) = DBLE(iPhase(1,iDCRT(lDCRT)))*Ccoor(1)
         TC(2) = DBLE(iPhase(2,iDCRT(lDCRT)))*Ccoor(2)
         TC(3) = DBLE(iPhase(3,iDCRT(lDCRT)))*Ccoor(3)


C Generate the quadrupole integral tables:
         iComp=6
         iOrdOp = 2
         nHer = (la + (lb+2) + 2 + 2) / 2
         Call MltPrm(Alpha,nAlpha,Beta,nBeta,Zeta,ZInv,rKappa,P,
     &             Array(ipTpp),nZeta,iComp,la,lb+2,A,RB,nHer,
     &             Array(ipArr),mArr,TC,iOrdOp)
         nHer = (la +  lb    + 2 + 2) / 2
         Call MltPrm(Alpha,nAlpha,Beta,nBeta,Zeta,ZInv,rKappa,P,
     &             Array(ipT  ),nZeta,iComp,la,lb  ,A,RB,nHer,
     &             Array(ipArr),mArr,TC,iOrdOp)
         if(lb.ge.2) then
           nHer = (la + (lb-2) + 2 + 2) / 2
           Call MltPrm(Alpha,nAlpha,Beta,nBeta,Zeta,ZInv,rKappa,P,
     &               Array(ipTmm),nZeta,iComp,la,lb-2,A,RB,nHer,
     &               Array(ipArr),mArr,TC,iOrdOp)
         end if
C Generate the dipole integral tables:
         iComp=3
         iOrdOp = 1
         nHer = (la + (lb+1) + 1 + 2) / 2
         Call MltPrm(Alpha,nAlpha,Beta,nBeta,Zeta,ZInv,rKappa,P,
     &             Array(ipTp ),nZeta,iComp,la,lb+1,A,RB,nHer,
     &             Array(ipArr),mArr,TC,iOrdOp)
         if(lb.ge.1) then
           nHer = (la + (lb-1) + 1 + 2) / 2
           Call MltPrm(Alpha,nAlpha,Beta,nBeta,Zeta,ZInv,rKappa,P,
     &               Array(ipTm ),nZeta,iComp,la,lb-1,A,RB,nHer,
     &               Array(ipArr),mArr,TC,iOrdOp)
         end if

         if(iprint.gt.49) write(6,*)' AMPInt calling AMPr.'
         Call AMPr(Array(ipB),nZeta,Array(ipRes),la,lb,Array(ipTpp),
     &             Array(ipTp),Array(ipT),Array(ipTm),Array(ipTmm))

C Symmetry adaption:
         if(iprint.gt.49) write(6,*)' AMPInt calling SymAdO'
         nOp = NrOpr(iDCRT(lDCRT),iOper,nIrrep)
         Call SymAdO(Array(ipRes),nZeta,la,lb,nComp,Final,nIC,
     &               nOp,lOper,iChO,One)
         if(iprint.gt.49) write(6,*)' Back to AMPInt.'
      End Do

      Call qExit('AMPInt')
      if(iprint.gt.49) write(6,*)' Leaving AMPInt.'
      Return
c Avoid unused argument warnings
      If (.False.) Then
         Call Unused_integer(nOrdOp)
         Call Unused_real(PtChrg)
         Call Unused_integer(nGrid)
         Call Unused_integer(iAddPot)
      End If
      End
