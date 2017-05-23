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
      SubRoutine CmbnM2(Rnxyz,nZeta,la,lb,Zeta,rKappa,Final,Alpha,Beta,
     &                  IfGrad,Fact,mVec)
************************************************************************
*                                                                      *
* Object: compute the gradient of the overlap matrix.                  *
*                                                                      *
* Called from: M2Grd                                                   *
*                                                                      *
* Calling    : QEnter                                                  *
*              DDot_   (ESSL)                                          *
*              QExit                                                   *
*                                                                      *
*     Author: Roland Lindh, Dept. of Theoretical Chemistry,            *
*             University of Lund, SWEDEN                               *
*             October '91.                                             *
************************************************************************
      Implicit Real*8 (A-H,O-Z)
#include "print.fh"
#include "real.fh"
      Real*8 Final(nZeta,(la+1)*(la+2)/2,(lb+1)*(lb+2)/2,6),
     &       Zeta(nZeta), rKappa(nZeta), Beta(nZeta),
     &       Rnxyz(nZeta,3,0:la+1,0:lb+1), Alpha(nZeta)
      Logical IfGrad(3,2)
*
*     Statement function for Cartesian index
*
      Ind(ixyz,ix,iz) = (ixyz-ix)*(ixyz-ix+1)/2 + iz + 1
*
      iRout = 134
      iPrint = nPrint(iRout)
*     Call qEnter('CmbnM2')
*     Call GetMem(' Enter CmbnM2','LIST','REAL',iDum,iDum)
*
*     ii = la*(la+1)*(la+2)/6
*     jj = lb*(lb+1)*(lb+2)/6
      exp32 = (-Three/Two)
      Do 25 iZeta = 1, nZeta
         rKappa(iZeta) = Fact* rKappa(iZeta) * Zeta(iZeta)**exp32
 25   Continue
      If (iPrint.ge.99) Then
         Call RecPrt(' In CmbnM2: Zeta  ',' ',Zeta  ,1,nZeta)
         Call RecPrt(' In CmbnM2: rKappa',' ',rKappa,1,nZeta)
         Call RecPrt(' In CmbnM2: Alpha ',' ',Alpha ,1,nZeta)
         Call RecPrt(' In CmbnM2: Beta  ',' ',Beta  ,1,nZeta)
         Call RecPrt(' In CmbnM2: Rnxyz ',' ',Rnxyz ,nZeta*3,
     &               (la+2)*(lb+2))
      End If
      Do 10 ixa = 0, la
         iyaMax=la-ixa
      Do 10 ixb = 0, lb
         iybMax=lb-ixb
         Do 20 iya = 0, iyaMax
            iza = la-ixa-iya
            ipa= Ind(la,ixa,iza)
         Do 20 iyb = 0, iybMax
            izb = lb-ixb-iyb
            ipb= Ind(lb,ixb,izb)
*
*           Combine overlap integrals
*
            tTwo = Two
            mVec = 0
            If (IfGrad(1,1)) Then
               mVec = mVec + 1
               If (ixa.gt.0) Then
                  xa = -DBLE(ixa)
                  Do 30 iZeta = 1, nZeta
                     Final(iZeta,ipa,ipb,mVec) = rKappa(iZeta)*
     &    (tTwo*Alpha(iZeta)*Rnxyz(iZeta,1,ixa+1,ixb) +
     &                    xa*Rnxyz(iZeta,1,ixa-1,ixb))*
     &                       Rnxyz(iZeta,2,iya,iyb)*
     &                       Rnxyz(iZeta,3,iza,izb)
 30               Continue
               Else
                  Do 31 iZeta = 1, nZeta
                     Final(iZeta,ipa,ipb,mVec) = rKappa(iZeta)*
     &     tTwo*Alpha(iZeta)*Rnxyz(iZeta,1,ixa+1,ixb)*
     &                       Rnxyz(iZeta,2,iya,iyb)*
     &                       Rnxyz(iZeta,3,iza,izb)
 31               Continue
               End If
            End If
            If (IfGrad(1,2)) Then
               mVec = mVec + 1
               If (ixb.gt.0) Then
                  xb = -DBLE(ixb)
                  Do 35 iZeta = 1, nZeta
                     Final(iZeta,ipa,ipb,mVec) = rKappa(iZeta)*
     &     (tTwo*Beta(iZeta)*Rnxyz(iZeta,1,ixa,ixb+1) +
     &                    xb*Rnxyz(iZeta,1,ixa,ixb-1))*
     &                       Rnxyz(iZeta,2,iya,iyb)*
     &                       Rnxyz(iZeta,3,iza,izb)
 35               Continue
               Else
                  Do 36 iZeta = 1, nZeta
                     Final(iZeta,ipa,ipb,mVec) = rKappa(iZeta)*
     &      tTwo*Beta(iZeta)*Rnxyz(iZeta,1,ixa,ixb+1)*
     &                       Rnxyz(iZeta,2,iya,iyb)*
     &                       Rnxyz(iZeta,3,iza,izb)
 36               Continue
               End If
            End If
            If (IfGrad(2,1)) Then
               mVec = mVec + 1
               If (iya.gt.0) Then
                  ya = -DBLE(iya)
                  Do 40 iZeta = 1, nZeta
                     Final(iZeta,ipa,ipb,mVec) = rKappa(iZeta)*
     &                       Rnxyz(iZeta,1,ixa,ixb)*
     &    (tTwo*Alpha(iZeta)*Rnxyz(iZeta,2,iya+1,iyb) +
     &                    ya*Rnxyz(iZeta,2,iya-1,iyb))*
     &                       Rnxyz(iZeta,3,iza,izb)
 40               Continue
               Else
                  Do 41 iZeta = 1, nZeta
                     Final(iZeta,ipa,ipb,mVec) = rKappa(iZeta)*
     &                       Rnxyz(iZeta,1,ixa,ixb)*
     &     tTwo*Alpha(iZeta)*Rnxyz(iZeta,2,iya+1,iyb)*
     &                       Rnxyz(iZeta,3,iza,izb)
 41               Continue
               End If
            End If
            If (IfGrad(2,2)) Then
               mVec = mVec + 1
               If (iyb.gt.0) Then
                  yb = -DBLE(iyb)
                  Do 45 iZeta = 1, nZeta
                     Final(iZeta,ipa,ipb,mVec) = rKappa(iZeta)*
     &                       Rnxyz(iZeta,1,ixa,ixb)*
     &     (tTwo*Beta(iZeta)*Rnxyz(iZeta,2,iya,iyb+1) +
     &                    yb*Rnxyz(iZeta,2,iya,iyb-1))*
     &                       Rnxyz(iZeta,3,iza,izb)
 45               Continue
               Else
                  Do 46 iZeta = 1, nZeta
                     Final(iZeta,ipa,ipb,mVec) = rKappa(iZeta)*
     &                       Rnxyz(iZeta,1,ixa,ixb)*
     &      tTwo*Beta(iZeta)*Rnxyz(iZeta,2,iya,iyb+1)*
     &                       Rnxyz(iZeta,3,iza,izb)
 46               Continue
               End If
            End If
            If (IfGrad(3,1)) Then
               mVec = mVec + 1
               If (iza.gt.0) Then
                  za = -DBLE(iza)
                  Do 50 iZeta = 1, nZeta
                     Final(iZeta,ipa,ipb,mVec) = rKappa(iZeta)*
     &                       Rnxyz(iZeta,1,ixa,ixb)*
     &                       Rnxyz(iZeta,2,iya,iyb)*
     &    (tTwo*Alpha(iZeta)*Rnxyz(iZeta,3,iza+1,izb) +
     &                    za*Rnxyz(iZeta,3,iza-1,izb))
 50               Continue
               Else
                  Do 51 iZeta = 1, nZeta
                     Final(iZeta,ipa,ipb,mVec) = rKappa(iZeta)*
     &                       Rnxyz(iZeta,1,ixa,ixb)*
     &                       Rnxyz(iZeta,2,iya,iyb)*
     &     tTwo*Alpha(iZeta)*Rnxyz(iZeta,3,iza+1,izb)
 51               Continue
               End If
            End If
            If (IfGrad(3,2)) Then
               mVec = mVec + 1
               If (izb.gt.0) Then
                  zb = -DBLE(izb)
                  Do 55 iZeta = 1, nZeta
                     Final(iZeta,ipa,ipb,mVec) = rKappa(iZeta)*
     &                       Rnxyz(iZeta,1,ixa,ixb)*
     &                       Rnxyz(iZeta,2,iya,iyb)*
     &     (tTwo*Beta(iZeta)*Rnxyz(iZeta,3,iza,izb+1) +
     &                    zb*Rnxyz(iZeta,3,iza,izb-1))
 55               Continue
               Else
                  Do 56 iZeta = 1, nZeta
                     Final(iZeta,ipa,ipb,mVec) = rKappa(iZeta)*
     &                       Rnxyz(iZeta,1,ixa,ixb)*
     &                       Rnxyz(iZeta,2,iya,iyb)*
     &      tTwo*Beta(iZeta)*Rnxyz(iZeta,3,iza,izb+1)
 56               Continue
               End If
            End If
*
 20      Continue
 10   Continue
*
*     Call GetMem(' Exit CmbnM2','LIST','REAL',iDum,iDum)
*     Call qExit('CmbnM2')
      Return
      End
