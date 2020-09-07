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
      Subroutine Darwin(Zeta,P,nZeta,A,Axyz,la,RB,Bxyz,lb, Final,
     &                  iStabM,nStabM,nComp,rKappa)
************************************************************************
*                                                                      *
* Object: to compoute the 1-electron Darwin contact term.              *
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
      use Basis_Info
      use Center_Info
      Implicit Real*8 (A-H,O-Z)
#include "itmax.fh"
#include "info.fh"
#include "constants.fh"
#include "WrkSpc.fh"
#include "print.fh"
#include "real.fh"
      Real*8 Final(nZeta,(la+1)*(la+2)/2,(lb+1)*(lb+2)/2,nComp),
     &       rKappa(nZeta),
     &       Axyz(nZeta,3,0:la), Bxyz(nZeta,3,0:lb),
     &       Zeta(nZeta), P(nZeta,3), A(3), RB(3), C(3), TC(3)
      Integer iStabM(0:nStabM-1), iDCRT(0:7)
*
*     Statement function for Cartesian index
*
      nElem(ixyz) = (ixyz+1)*(ixyz+2)/2
      Ind(ixyz,ix,iz) = (ixyz-ix)*(ixyz-ix+1)/2 + iz + 1
*
      iRout = 170
      iPrint = nPrint(iRout)
      Call qEnter('Darwin ')
      If (iPrint.ge.99) Then
         Call RecPrt(' In Darwin: rKappa',' ',rKappa,nZeta,1)
         Call RecPrt(' In Darwin: Zeta',' ',Zeta,nZeta,1)
         Call RecPrt(' In Darwin: P',' ',P,nZeta,3)
      End If
*
      call dcopy_(nZeta*nElem(la)*nElem(lb)*nComp,[Zero],0,Final,1)
*
      kdc = 0
      Do 500 kCnttp = 1, nCnttp
         Do 501 kCnt = 1, dbsc(kCnttp)%nCntr
            C(1:3) = dbsc(kCnttp)%Coor(1:3,kCnt)
*
            Call DCR(LmbdT,iOper,nIrrep,iStabM,nStabM,
     &               dc(kdc+kCnt)%iStab,dc(kdc+kCnt)%nStab,iDCRT,nDCRT)
            Fact = DBLE(nStabM) / DBLE(LmbdT)
*
            Do 502 lDCRT = 0, nDCRT-1
               Call OA(iDCRT(lDCRT),C,TC)
*
*--------------Compute the value of the angular components associated
*              to the basis functions centered on the first center.
*
               call dcopy_(nZeta*3,[One],0,Axyz(1,1,0),1)
               If (la.eq.0) Go To 60
*
               Do 20 iCar = 1, 3
*
                  Do 30 iZeta = 1, nZeta
                     Axyz(iZeta,iCar,1) =TC(iCar) -A(iCar)
 30               Continue
*
                  Do 40 ia = 2, la
                     Do 50 iZeta = 1, nZeta
                        Axyz(iZeta,iCar,ia) = Axyz(iZeta,iCar,1) *
     &                                        Axyz(iZeta,iCar,ia-1)
 50                  Continue
 40               Continue
*
 20            Continue
 60            Continue
*
*--------------Compute the value of the angular components associated to
*              the basis functions centered on the second center.
*
               call dcopy_(nZeta*3,[One],0,Bxyz(1,1,0),1)
*
*--------------Modify z-component to carry the charge and the exponetial
*              contribution.
*
               Do 210 iZeta = 1, nZeta
                  Bxyz(iZeta,3,0) = dbsc(kCnttp)%Charge *
     &                Exp ( -Zeta(iZeta) * ( (TC(1)-P(iZeta,1))**2 +
     &                                       (TC(2)-P(iZeta,2))**2 +
     &                                       (TC(3)-P(iZeta,3))**2 ))
 210           Continue
               If (lb.eq.0) Go To 61
*
               Do 21 iCar = 1, 3
*
                  Do 31 iZeta = 1, nZeta
                     Bxyz(iZeta,iCar,1) =TC(iCar) -RB(iCar)
 31               Continue
*
                  Do 41 ib = 2, lb
                     Do 51 iZeta = 1, nZeta
                        Bxyz(iZeta,iCar,ib) = Bxyz(iZeta,iCar,1) *
     &                                        Bxyz(iZeta,iCar,ib-1)
 51                  Continue
 41               Continue
 21            Continue
*
*--------------Modify z-components with the exponential contribution
*
               Do 32 ib = 1, lb
                  Do 42 iZeta = 1, nZeta
                     Bxyz(iZeta,3,ib) = Bxyz(iZeta,3,ib) *
     &                                  Bxyz(iZeta,3,0)
 42               Continue
 32            Continue
*
 61            Continue
*
*--------------Combine contributions from the varoius angular
*              components.
*
               Do 100 ixa = la, 0, -1
               Do 101 ixb = lb, 0, -1
                  Do 110 iya = la-ixa, 0, -1
                     iza = la-ixa-iya
                     ipa = Ind(la,ixa,iza)
                  Do 111 iyb = lb-ixb, 0, -1
                     izb = lb-ixb-iyb
                     ipb = Ind(lb,ixb,izb)
                     Do 130 iZeta = 1, nZeta
                        Final(iZeta,ipa,ipb,1) =
     &                     Final(iZeta,ipa,ipb,1) + Fact *
     &                                Axyz(iZeta,1,ixa) *
     &                                Axyz(iZeta,2,iya) *
     &                                Axyz(iZeta,3,iza) *
     &                                Bxyz(iZeta,1,ixb) *
     &                                Bxyz(iZeta,2,iyb) *
     &                                Bxyz(iZeta,3,izb)
 130                 Continue
 111              Continue
 110              Continue
 101           Continue
 100           Continue

 502        Continue
 501     Continue
         kdc = kdc + dbsc(kCnttp)%nCntr
 500  Continue
*
*     Factor from operator (pi/(2*c**2), c=137.036 au)
*
*     Factor = Pi * One2C2
      Factor = Pi / (Two* CONST_C_IN_AU_ **2)
      Do 140 ipa = 1, nElem(la)
         Do 141 ipb = 1, nELem(lb)
            Do 142 iZeta = 1, nZeta
               Final(iZeta,ipa,ipb,1) = rKappa(iZeta) *
     &          Factor * Final(iZeta,ipa,ipb,1)
 142        Continue
 141     Continue
 140  Continue
*
*     Call GetMem(' Exit Darwin ','LIST','REAL',iDum,iDum)
      Call qExit('Darwin ')
      Return
      End
