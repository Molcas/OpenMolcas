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
* Copyright (C) 1995, Roland Lindh                                     *
************************************************************************
      SubRoutine CmbnW1(Welp0,Welm0,Wel0p,Wel0m,
     &                  nZeta,la,lb,Zeta,rKappa,Final,Alpha,Beta,
     &                  Grad,nGrad,DAO,IfGrad,IndGrd,iStab,jStab,kOp)
************************************************************************
*                                                                      *
* Object: compute the gradient of the Spherical Well integrals         *
*                                                                      *
* Called from: WelGrd                                                  *
*                                                                      *
* Calling    : QEnter                                                  *
*              DDot_   (ESSL)                                          *
*              QExit                                                   *
*                                                                      *
*     Author: Roland Lindh, Dept. of Theoretical Chemistry,            *
*             University of Lund, SWEDEN                               *
*             May '95.                                                 *
************************************************************************
      use Symmetry_Info, only: nIrrep,iChBas
      Implicit Real*8 (A-H,O-Z)
#include "print.fh"
#include "real.fh"
      Real*8 Final(nZeta,(la+1)*(la+2)/2,(lb+1)*(lb+2)/2,6),
     &       Zeta(nZeta), rKappa(nZeta), Beta(nZeta),
     &       Welp0(nZeta,(la+2)*(la+3)/2,(lb+1)*(lb+2)/2),
     &       Welm0(nZeta, la   *(la+1)/2,(lb+1)*(lb+2)/2),
     &       Wel0p(nZeta,(la+1)*(la+2)/2,(lb+2)*(lb+3)/2),
     &       Wel0m(nZeta,(la+1)*(la+2)/2, lb   *(lb+1)/2),
     &       Alpha(nZeta), Grad(nGrad),
     &       DAO(nZeta,(la+1)*(la+2)/2,(lb+1)*(lb+2)/2)
      Logical IfGrad(3,2)
      Integer IndGrd(3,2), kOp(2)
*
*     Statement function for Cartesian index
*
*                     ix*0   added to avoid compiler warning
      Ind(ix,iy,iz) = ix*0 + (iy+iz)*(iy+iz+1)/2 + iz + 1
*
      iRout = 134
      iPrint = nPrint(iRout)
      iQ = 1
*     Call qEnter('CmbnW1')
*
      If (iPrint.ge.99) Then
         Call RecPrt(' In CmbnW1: Zeta  ',' ',Zeta  ,1,nZeta)
         Call RecPrt(' In CmbnW1: rKappa',' ',rKappa,1,nZeta)
         Call RecPrt(' In CmbnW1: Alpha ',' ',Alpha ,1,nZeta)
         Call RecPrt(' In CmbnW1: Beta  ',' ',Beta  ,1,nZeta)
         np0=(la+2)*(la+3)/2*(lb+1)*(lb+2)/2
         nm0= la   *(la+1)/2*(lb+1)*(lb+2)/2
         n0p=(la+1)*(la+2)/2*(lb+2)*(lb+3)/2
         n0m=(la+1)*(la+2)/2* lb   *(lb+1)/2
         Call RecPrt(' In CmbnW1: Welp0',' ',Welp0,nZeta,np0)
         If (la.ge.1)
     &      Call RecPrt(' In CmbnW1: Welm0',' ',Welm0,nZeta,nm0)
         Call RecPrt(' In CmbnW1: Wel0p',' ',Wel0p,nZeta,n0p)
         If (lb.ge.1)
     &      Call RecPrt(' In CmbnW1: Wel0m',' ',Wel0m,nZeta,n0m)
      End If
*
      tTwo = Two
      nDAO = nZeta * (la+1)*(la+2)/2 * (lb+1)*(lb+2)/2
      Do 10 ixa = 0, la
         iyaMax=la-ixa
      Do 11 ixb = 0, lb
         iybMax=lb-ixb
         Do 20 iya = 0, iyaMax
            iza = la-ixa-iya
            ipa= Ind(ixa,iya,iza)
         Do 21 iyb = 0, iybMax
            izb = lb-ixb-iyb
            ipb= Ind(ixb,iyb,izb)
*
*           Combine Spherical Well integrals
*
            If (IfGrad(1,1)) Then
               ipp0=Ind(ixa+1,iya,iza)
               If (ixa.gt.0) Then
                  xa = -ixa
                  ipm0=Ind(ixa-1,iya,iza)
                  Do iZeta = 1, nZeta
                     Final(iZeta,ipa,ipb,1) =
     &                   (tTwo*Alpha(iZeta)*Welp0(iZeta,ipp0,ipb) +
     &                                   xa*Welm0(iZeta,ipm0,ipb) )
                  End Do
               Else
                  Do iZeta = 1, nZeta
                     Final(iZeta,ipa,ipb,1) =
     &                   (tTwo*Alpha(iZeta)*Welp0(iZeta,ipp0,ipb) )
                  End Do
               End If
            End If
            If (IfGrad(1,2)) Then
               ip0p=Ind(ixb+1,iyb,izb)
               If (ixb.gt.0) Then
                  xb = -ixb
                  ip0m=Ind(ixb-1,iyb,izb)
                  Do iZeta = 1, nZeta
                     Final(iZeta,ipa,ipb,4) =
     &                   (tTwo* Beta(iZeta)*Wel0p(iZeta,ipa,ip0p) +
     &                                   xb*Wel0m(iZeta,ipa,ip0m) )
                  End Do
               Else
                  Do iZeta = 1, nZeta
                     Final(iZeta,ipa,ipb,4) =
     &                   (tTwo* Beta(iZeta)*Wel0p(iZeta,ipa,ip0p) )
                  End Do
               End If
            End If
            If (IfGrad(2,1)) Then
               ipp0=Ind(ixa,iya+1,iza)
               If (iya.gt.0) Then
                  ya = -iya
                  ipm0=Ind(ixa,iya-1,iza)
                  Do iZeta = 1, nZeta
                     Final(iZeta,ipa,ipb,2) =
     &                   (tTwo*Alpha(iZeta)*Welp0(iZeta,ipp0,ipb) +
     &                                   ya*Welm0(iZeta,ipm0,ipb) )
                  End Do
               Else
                  Do iZeta = 1, nZeta
                     Final(iZeta,ipa,ipb,2) =
     &                   (tTwo*Alpha(iZeta)*Welp0(iZeta,ipp0,ipb) )
                  End Do
               End If
            End If
            If (IfGrad(2,2)) Then
               ip0p=Ind(ixb,iyb+1,izb)
               If (iyb.gt.0) Then
                  yb = -iyb
                  ip0m=Ind(ixb,iyb-1,izb)
                  Do iZeta = 1, nZeta
                     Final(iZeta,ipa,ipb,5) =
     &                   (tTwo* Beta(iZeta)*Wel0p(iZeta,ipa,ip0p) +
     &                                   yb*Wel0m(iZeta,ipa,ip0m) )
                  End Do
               Else
                  Do iZeta = 1, nZeta
                     Final(iZeta,ipa,ipb,5) =
     &                   (tTwo* Beta(iZeta)*Wel0p(iZeta,ipa,ip0p) )
                  End Do
               End If
            End If
            If (IfGrad(3,1)) Then
               ipp0=Ind(ixa,iya,iza+1)
               If (iza.gt.0) Then
                  za = -iza
                  ipm0=Ind(ixa,iya,iza-1)
                  Do iZeta = 1, nZeta
                     Final(iZeta,ipa,ipb,3) =
     &                   (tTwo*Alpha(iZeta)*Welp0(iZeta,ipp0,ipb) +
     &                                   za*Welm0(iZeta,ipm0,ipb) )
                  End Do
               Else
                  Do iZeta = 1, nZeta
                     Final(iZeta,ipa,ipb,3) =
     &                   (tTwo*Alpha(iZeta)*Welp0(iZeta,ipp0,ipb) )
                  End Do
               End If
            End If
            If (IfGrad(3,2)) Then
               ip0p=Ind(ixb,iyb,izb+1)
               If (izb.gt.0) Then
                  zb = -izb
                  ip0m=Ind(ixb,iyb,izb-1)
                  Do iZeta = 1, nZeta
                     Final(iZeta,ipa,ipb,6) =
     &                   (tTwo* Beta(iZeta)*Wel0p(iZeta,ipa,ip0p) +
     &                                   zb*Wel0m(iZeta,ipa,ip0m) )
                  End Do
               Else
                  Do iZeta = 1, nZeta
                     Final(iZeta,ipa,ipb,6) =
     &                   (tTwo* Beta(iZeta)*Wel0p(iZeta,ipa,ip0p) )
                  End Do
               End If
            End If
*
 21      Continue
 20      Continue
 11   Continue
 10   Continue
*
*     Trace the gradient integrals
*
      If (iPrint.ge.99) Then
         Call RecPrt(' W(1)',' ',Final,nDAO,6)
         Call RecPrt('   D ',' ',DAO,nDAO,1)
      End If
      Do 100 iCn = 1, 2
         Do 110 iCar = 1, 3
            If (IndGrd(iCar,iCn).ne.0) Then
*              Accumulate contribution to the gradient
               iGrad = Abs(IndGrd(iCar,iCn))
               If (iCn.eq.1) Then
                   i1 = iCar
                   i2 = iCar + 3
                   ps = DBLE( iPrmt( kOp(1), iChBas(1+iCar) ) )
                   Fact = DBLE(iStab)/DBLE(nIrrep)
               Else
                   i1 = iCar + 3
                   i2 = iCar
                   ps = DBLE( iPrmt( kOp(2), iChBas(1+iCar) ) )
                   Fact = ps * DBLE(jStab)/DBLE(nIrrep)
               End if
*              Write (*,*) ' Fact=', Fact, ps
               If (IndGrd(iCar,iCn).lt.0) Then
*                 Gradient via translational invariance.
                  Grad(iGrad) = Grad(iGrad) - Fact*
     &                DDot_(nDAO,DAO,1,Final(1,1,1,i2),1)
               Else
                  Grad(iGrad) = Grad(iGrad) + Fact*
     &                DDot_(nDAO,DAO,1,Final(1,1,1,i1),1)
               End If
            End If
 110     Continue
 100  Continue
*
*     Call qExit('CmbnW1')
      Return
      End
