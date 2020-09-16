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
      SubRoutine CmbnS1(Rnxyz,nZeta,la,lb,Zeta,rKappa,Final,Alpha,Beta,
     &                  Grad,nGrad,DAO,IfGrad,IndGrd,iStab,jStab,kOp)
************************************************************************
*                                                                      *
* Object: compute the gradient of the overlap matrix.                  *
*                                                                      *
* Called from: OvrGrd                                                  *
*                                                                      *
* Calling    : QEnter                                                  *
*              DDot_   (ESSL)                                          *
*              QExit                                                   *
*                                                                      *
*     Author: Roland Lindh, Dept. of Theoretical Chemistry,            *
*             University of Lund, SWEDEN                               *
*             October '91.                                             *
************************************************************************
      use Symmetry_Info, only: iChBas, nIrrep
      Implicit Real*8 (A-H,O-Z)
#include "print.fh"
#include "real.fh"
#include "nac.fh"
      Real*8 Final(nZeta,(la+1)*(la+2)/2,(lb+1)*(lb+2)/2,6),
     &       Zeta(nZeta), rKappa(nZeta), Beta(nZeta),
     &       Rnxyz(nZeta,3,0:la+1,0:lb+1), Alpha(nZeta), Grad(nGrad),
     &       DAO(nZeta,(la+1)*(la+2)/2,(lb+1)*(lb+2)/2)
      Logical IfGrad(3,2)
      Integer IndGrd(3,2), kOp(2)
*
*     Statement function for Cartesian index
*
      Ind(ixyz,ix,iz) = (ixyz-ix)*(ixyz-ix+1)/2 + iz + 1
*
      iRout = 134
      iPrint = nPrint(iRout)
      iQ = 1
*     Call qEnter('CmbnS1')
*
*     ii = la*(la+1)*(la+2)/6
*     jj = lb*(lb+1)*(lb+2)/6
      exp32 = -Three/Two
      Do 25 iZeta = 1, nZeta
         rKappa(iZeta) = rKappa(iZeta) * Zeta(iZeta)**exp32
 25   Continue
      If (iPrint.ge.99) Then
         Call RecPrt(' In CmbnS1: Zeta  ',' ',Zeta  ,1,nZeta)
         Call RecPrt(' In CmbnS1: rKappa',' ',rKappa,1,nZeta)
         Call RecPrt(' In CmbnS1: Alpha ',' ',Alpha ,1,nZeta)
         Call RecPrt(' In CmbnS1: Beta  ',' ',Beta  ,1,nZeta)
      End If
      Do 10 ixa = 0, la
         iyaMax=la-ixa
      Do 11 ixb = 0, lb
         iybMax=lb-ixb
         Do 20 iya = 0, iyaMax
            iza = la-ixa-iya
            ipa= Ind(la,ixa,iza)
*           iChBs = iChBas(ii+ipa)
*           pa    = DBLE(iPrmt(kOp(1),iChBs))
         Do 21 iyb = 0, iybMax
            izb = lb-ixb-iyb
            ipb= Ind(lb,ixb,izb)
*           jChBs = iChBas(jj+ipb)
*           papb  = DBLE(iPrmt(kOp(2),jChBs)) * pa
*
*           Combine overlap integrals
*
            tTwo = Two
*           write (*,*) ' papb=', papb
            If (IfGrad(1,1)) Then
               If (ixa.gt.0) Then
                  xa = Dble(-ixa)
                  Do 30 iZeta = 1, nZeta
*                    Final(iZeta,ipa,ipb,1) = papb * rKappa(iZeta)*
                     Final(iZeta,ipa,ipb,1) = rKappa(iZeta)*
     &    (tTwo*Alpha(iZeta)*Rnxyz(iZeta,1,ixa+1,ixb) +
     &                    xa*Rnxyz(iZeta,1,ixa-1,ixb))*
     &                       Rnxyz(iZeta,2,iya,iyb)*
     &                       Rnxyz(iZeta,3,iza,izb)
 30               Continue
               Else
                  Do 31 iZeta = 1, nZeta
*                    Final(iZeta,ipa,ipb,1) = papb * rKappa(iZeta)*
                     Final(iZeta,ipa,ipb,1) = rKappa(iZeta)*
     &     tTwo*Alpha(iZeta)*Rnxyz(iZeta,1,ixa+1,ixb)*
     &                       Rnxyz(iZeta,2,iya,iyb)*
     &                       Rnxyz(iZeta,3,iza,izb)
 31               Continue
               End If
            End If
            If (IfGrad(1,2)) Then
               If (ixb.gt.0) Then
                  xb = Dble(-ixb)
                  Do 35 iZeta = 1, nZeta
*                    Final(iZeta,ipa,ipb,4) = papb * rKappa(iZeta)*
                     Final(iZeta,ipa,ipb,4) = rKappa(iZeta)*
     &     (tTwo*Beta(iZeta)*Rnxyz(iZeta,1,ixa,ixb+1) +
     &                    xb*Rnxyz(iZeta,1,ixa,ixb-1))*
     &                       Rnxyz(iZeta,2,iya,iyb)*
     &                       Rnxyz(iZeta,3,iza,izb)
 35               Continue
               Else
                  Do 36 iZeta = 1, nZeta
*                    Final(iZeta,ipa,ipb,4) = papb * rKappa(iZeta)*
                     Final(iZeta,ipa,ipb,4) = rKappa(iZeta)*
     &      tTwo*Beta(iZeta)*Rnxyz(iZeta,1,ixa,ixb+1)*
     &                       Rnxyz(iZeta,2,iya,iyb)*
     &                       Rnxyz(iZeta,3,iza,izb)
 36               Continue
               End If
            End If
            If (IfGrad(2,1)) Then
               If (iya.gt.0) Then
                  ya = Dble(-iya)
                  Do 40 iZeta = 1, nZeta
*                    Final(iZeta,ipa,ipb,2) = papb * rKappa(iZeta)*
                     Final(iZeta,ipa,ipb,2) = rKappa(iZeta)*
     &                       Rnxyz(iZeta,1,ixa,ixb)*
     &    (tTwo*Alpha(iZeta)*Rnxyz(iZeta,2,iya+1,iyb) +
     &                    ya*Rnxyz(iZeta,2,iya-1,iyb))*
     &                       Rnxyz(iZeta,3,iza,izb)
 40               Continue
               Else
                  Do 41 iZeta = 1, nZeta
*                    Final(iZeta,ipa,ipb,2) = papb * rKappa(iZeta)*
                     Final(iZeta,ipa,ipb,2) = rKappa(iZeta)*
     &                       Rnxyz(iZeta,1,ixa,ixb)*
     &     tTwo*Alpha(iZeta)*Rnxyz(iZeta,2,iya+1,iyb)*
     &                       Rnxyz(iZeta,3,iza,izb)
 41               Continue
               End If
            End If
            If (IfGrad(2,2)) Then
               If (iyb.gt.0) Then
                  yb = Dble(-iyb)
                  Do 45 iZeta = 1, nZeta
*                    Final(iZeta,ipa,ipb,5) = papb * rKappa(iZeta)*
                     Final(iZeta,ipa,ipb,5) = rKappa(iZeta)*
     &                       Rnxyz(iZeta,1,ixa,ixb)*
     &     (tTwo*Beta(iZeta)*Rnxyz(iZeta,2,iya,iyb+1) +
     &                    yb*Rnxyz(iZeta,2,iya,iyb-1))*
     &                       Rnxyz(iZeta,3,iza,izb)
 45               Continue
               Else
                  Do 46 iZeta = 1, nZeta
*                    Final(iZeta,ipa,ipb,5) = papb * rKappa(iZeta)*
                     Final(iZeta,ipa,ipb,5) = rKappa(iZeta)*
     &                       Rnxyz(iZeta,1,ixa,ixb)*
     &      tTwo*Beta(iZeta)*Rnxyz(iZeta,2,iya,iyb+1)*
     &                       Rnxyz(iZeta,3,iza,izb)
 46               Continue
               End If
            End If
            If (IfGrad(3,1)) Then
               If (iza.gt.0) Then
                  za = Dble(-iza)
                  Do 50 iZeta = 1, nZeta
*                    Final(iZeta,ipa,ipb,3) = papb * rKappa(iZeta)*
                     Final(iZeta,ipa,ipb,3) = rKappa(iZeta)*
     &                       Rnxyz(iZeta,1,ixa,ixb)*
     &                       Rnxyz(iZeta,2,iya,iyb)*
     &    (tTwo*Alpha(iZeta)*Rnxyz(iZeta,3,iza+1,izb) +
     &                    za*Rnxyz(iZeta,3,iza-1,izb))
 50               Continue
               Else
                  Do 51 iZeta = 1, nZeta
*                    Final(iZeta,ipa,ipb,3) = papb * rKappa(iZeta)*
                     Final(iZeta,ipa,ipb,3) = rKappa(iZeta)*
     &                       Rnxyz(iZeta,1,ixa,ixb)*
     &                       Rnxyz(iZeta,2,iya,iyb)*
     &     tTwo*Alpha(iZeta)*Rnxyz(iZeta,3,iza+1,izb)
 51               Continue
               End If
            End If
            If (IfGrad(3,2)) Then
               If (izb.gt.0) Then
                  zb = Dble(-izb)
                  Do 55 iZeta = 1, nZeta
*                    Final(iZeta,ipa,ipb,6) = papb * rKappa(iZeta)*
                     Final(iZeta,ipa,ipb,6) = rKappa(iZeta)*
     &                       Rnxyz(iZeta,1,ixa,ixb)*
     &                       Rnxyz(iZeta,2,iya,iyb)*
     &     (tTwo*Beta(iZeta)*Rnxyz(iZeta,3,iza,izb+1) +
     &                    zb*Rnxyz(iZeta,3,iza,izb-1))
 55               Continue
               Else
                  Do 56 iZeta = 1, nZeta
*                    Final(iZeta,ipa,ipb,6) = papb * rKappa(iZeta)*
                     Final(iZeta,ipa,ipb,6) = rKappa(iZeta)*
     &                       Rnxyz(iZeta,1,ixa,ixb)*
     &                       Rnxyz(iZeta,2,iya,iyb)*
     &      tTwo*Beta(iZeta)*Rnxyz(iZeta,3,iza,izb+1)
 56               Continue
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
      nDAO = nZeta * (la+1)*(la+2)/2 * (lb+1)*(lb+2)/2
      If (iPrint.ge.99) Then
         Call RecPrt(' S(1)',' ',Final,nDAO,6)
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
                  IF (isCSF) Then
*                    The CSF NADC component is not translation invariant,
*                    in fact, it is based on an antisymmetric matrix, so a + sign is needed here
                     Grad(iGrad) = Grad(iGrad) + Fact*
     &                   DDot_(nDAO,DAO,1,Final(1,1,1,i2),1)
                  Else
*                    Gradient via translational invariance.
                     Grad(iGrad) = Grad(iGrad) - Fact*
     &                   DDot_(nDAO,DAO,1,Final(1,1,1,i2),1)
                  End if
               Else
                  Grad(iGrad) = Grad(iGrad) + Fact*
     &                DDot_(nDAO,DAO,1,Final(1,1,1,i1),1)
               End If
            End If
 110     Continue
 100  Continue
*
*     Call qExit('CmbnS1')
      Return
      End
