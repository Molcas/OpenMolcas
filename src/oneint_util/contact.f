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
      Subroutine Contact(Zeta,P,nZeta,A,Axyz,la,RB,Bxyz,lb, Ccoor,
     &                   lOper,iChO,nIC,Array,nArr,Final,
     &                   iStabM,nStabM,nComp,rKappa)
************************************************************************
*                                                                      *
* Object: to compoute the 1-electron contact term.                     *
*                                                                      *
* Called from: D1Int                                                   *
*                                                                      *
* Calling    : QEnter                                                  *
*              RecPrt                                                  *
*              GetMem                                                  *
*              QExit                                                   *
*                                                                      *
*     Author: Roland Lindh, Dept. Of Theoretical Chemistry,            *
*             University of Lund, Sweden, February '91                 *
************************************************************************
      Implicit Real*8 (A-H,O-Z)
#include "itmax.fh"
#include "info.fh"
#include "WrkSpc.fh"
#include "print.fh"
#include "real.fh"
      Real*8 Final(nZeta,(la+1)*(la+2)/2,(lb+1)*(lb+2)/2,nIC),
     &       rKappa(nZeta), Ccoor(3),
     &       Array(nZeta,(la+1)*(la+2)/2,(lb+1)*(lb+2)/2),
     &       Axyz(nZeta,3,0:la), Bxyz(nZeta,3,0:lb),
     &       Zeta(nZeta), P(nZeta,3), A(3), RB(3), TC(3)
      Integer iStabM(0:nStabM-1), iStabO(0:7), iDCRT(0:7), lOper(nComp),
     &       iChO(nComp)
*                                                                      *
************************************************************************
*                                                                      *
*     Statement function for Cartesian index
*
      nElem(ixyz) = (ixyz+1)*(ixyz+2)/2
      Ind(ixyz,ix,iz) = (ixyz-ix)*(ixyz-ix+1)/2 + iz + 1
*                                                                      *
************************************************************************
*                                                                      *
      iRout = 170
      iPrint = nPrint(iRout)
      Call qEnter('Contact ')
      If (iPrint.ge.99) Then
         Call RecPrt(' In Contact: rKappa',' ',rKappa,nZeta,1)
         Call RecPrt(' In Contact: Zeta',' ',Zeta,nZeta,1)
         Call RecPrt(' In Contact: P',' ',P,nZeta,3)
      End If
*                                                                      *
************************************************************************
*                                                                      *
      llOper = lOper(1)
      Do iComp = 2, nComp
         llOper = iOr(llOper,lOper(iComp))
      End Do
      Call SOS(iStabO,nStabO,llOper)
      Call DCR(LmbdT,iOper,nIrrep,iStabM,nStabM,iStabO,nStabO,
     &         iDCRT,nDCRT)
*
      Do lDCRT = 0, nDCRT-1
         TC(1) = DBLE(iPhase(1,iDCRT(lDCRT)))*Ccoor(1)
         TC(2) = DBLE(iPhase(2,iDCRT(lDCRT)))*Ccoor(2)
         TC(3) = DBLE(iPhase(3,iDCRT(lDCRT)))*Ccoor(3)
*
         call dcopy_(nZeta*nElem(la)*nElem(lb),Zero,0,Array,1)
*
*--------Compute the value of the angular components associated
*        to the basis functions centered on the first center.
*
         call dcopy_(nZeta*3,One,0,Axyz(1,1,0),1)
         If (la.eq.0) Go To 60
*
         Do iCar = 1, 3
*
            Do iZeta = 1, nZeta
               Axyz(iZeta,iCar,1) =TC(iCar) -A(iCar)
            End Do
*
            Do ia = 2, la
               Do iZeta = 1, nZeta
                  Axyz(iZeta,iCar,ia) = Axyz(iZeta,iCar,1) *
     &                                  Axyz(iZeta,iCar,ia-1)
               End Do
            End Do
*
         End Do
 60      Continue
*
*--------Compute the value of the angular components associated to
*        the basis functions centered on the second center.
*
         call dcopy_(nZeta*3,One,0,Bxyz(1,1,0),1)
*
*--------Modify z-component to carry the the exponetial
*        contribution.
*
         Do iZeta = 1, nZeta
            Bxyz(iZeta,3,0) =
     &          Exp ( -Zeta(iZeta) * ( (TC(1)-P(iZeta,1))**2 +
     &                                 (TC(2)-P(iZeta,2))**2 +
     &                                 (TC(3)-P(iZeta,3))**2 ))
         End Do
         If (lb.eq.0) Go To 61
*
         Do iCar = 1, 3
*
            Do iZeta = 1, nZeta
               Bxyz(iZeta,iCar,1) =TC(iCar) -RB(iCar)
            End Do
*
            Do ib = 2, lb
               Do iZeta = 1, nZeta
                  Bxyz(iZeta,iCar,ib) = Bxyz(iZeta,iCar,1) *
     &                                  Bxyz(iZeta,iCar,ib-1)
               End Do
            End Do
         End Do
*
*--------Modify z-components with the exponential contribution
*
         Do ib = 1, lb
            Do iZeta = 1, nZeta
               Bxyz(iZeta,3,ib) = Bxyz(iZeta,3,ib) *
     &                            Bxyz(iZeta,3,0)
            End Do
         End Do
*
 61      Continue
*
*--------Combine contributions from the varoius angular
*        components.
*
         Do ixa = la, 0, -1
            Do ixb = lb, 0, -1
               Do iya = la-ixa, 0, -1
                  iza = la-ixa-iya
                  ipa = Ind(la,ixa,iza)
                  Do iyb = lb-ixb, 0, -1
                     izb = lb-ixb-iyb
                     ipb = Ind(lb,ixb,izb)
                     Do iZeta = 1, nZeta
                        Array(iZeta,ipa,ipb) =
     &                     Array(iZeta,ipa,ipb) +
     &                                rKappa(iZeta) *
     &                                Axyz(iZeta,1,ixa) *
     &                                Axyz(iZeta,2,iya) *
     &                                Axyz(iZeta,3,iza) *
     &                                Bxyz(iZeta,1,ixb) *
     &                                Bxyz(iZeta,2,iyb) *
     &                                Bxyz(iZeta,3,izb)
                     End Do
                  End Do
               End Do
            End Do
         End Do
*
*------- Accumulate contributions
*
         nOp = NrOpr(iDCRT(lDCRT),iOper,nIrrep)
         Call SymAdO(Array,nZeta,la,lb,nComp,Final,nIC,
     &               nOp         ,lOper,iChO,One)
*
      End Do
*
*     Call GetMem(' Exit Contact ','LIST','REAL',iDum,iDum)
      Call qExit('Contact ')
      Return
c Avoid unused argument warnings
      If (.False.) Call Unused_integer(nArr)
      End
