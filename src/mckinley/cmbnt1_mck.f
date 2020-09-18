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
*               1991, Anders Bernhardsson                              *
************************************************************************
      SubRoutine CmbnT1_mck(Rnxyz,nZeta,la,lb,Zeta,rKappa,
     &                  Final,Txyz,Alpha,Beta,IfGrad)
************************************************************************
*                                                                      *
* Object:                                                              *
*                                                                      *
* Called from: KnEGrd                                                  *
*                                                                      *
* Calling    : QEnter                                                  *
*              DDot_   (ESSL)                                          *
*              QExit                                                   *
*                                                                      *
*     Author: Roland Lindh, Dept. of Theoretical Chemistry,            *
*             University of Lund, SWEDEN                               *
*             October '91                                              *
*             Anders Bernhardsson, Dept. of Theoretical Chemistry,     *
*             University of Lund, SWEDEN                               *
*             October '91                                              *
************************************************************************
      Implicit Real*8 (A-H,O-Z)
c#include "print.fh"
#include "real.fh"
      Real*8 Final(nZeta,(la+1)*(la+2)/2,(lb+1)*(lb+2)/2,1),
     &       Zeta(nZeta), rKappa(nZeta), Alpha(nZeta), Beta(nZeta),
     &       Rnxyz(nZeta,3,0:la+2,0:lb+2),
     &       Txyz(nZeta,3,0:la+1,0:lb+1)
      Logical IfGrad(3,2)
*
*     Statement function for Cartesian index
*
      Ind(ixyz,ix,iz) = (ixyz-ix)*(ixyz-ix+1)/2 + iz + 1
*
c     iRout = 134
c     iPrint = nPrint(iRout)
c     Call QEnter('CmbnT1')
*
*     ii = la*(la+1)*(la+2)/6
*     jj = lb*(lb+1)*(lb+2)/6
      exp32 = -Three/Two
      Do 25 iZeta = 1, nZeta
         rKappa(iZeta) = rKappa(iZeta) * Zeta(iZeta)**exp32
 25   Continue
      Do 10 ixa = 0, la
         iyaMax=la-ixa
      Do 11 ixb = 0, lb
         iybMax=lb-ixb
         Do 20 iya = 0, iyaMax
            iza = la-ixa-iya
            ipa= Ind(la,ixa,iza)
         Do 21 iyb = 0, iybMax
            izb = lb-ixb-iyb
            ipb= Ind(lb,ixb,izb)
*
*           Combine integrals
*
            tTwo = Two
            If (IfGrad(1,1)) Then
               If (ixa.gt.0) Then
                  xa = Dble(-ixa)
                  Do 30 iZeta = 1, nZeta
                     Final(iZeta,ipa,ipb,1) =
     &                rKappa(iZeta) * (
     &          (tTwo*Txyz(iZeta,1,ixa+1,ixb)*Alpha(iZeta) +
     &             xa*Txyz(iZeta,1,ixa-1,ixb) )*
     &               Rnxyz(iZeta,2,iya,iyb)*
     &               Rnxyz(iZeta,3,iza,izb) +
     &         (tTwo*Rnxyz(iZeta,1,ixa+1,ixb)*Alpha(iZeta) +
     &            xa*Rnxyz(iZeta,1,ixa-1,ixb) )*
     &               Txyz(iZeta,2,iya,iyb)*
     &               Rnxyz(iZeta,3,iza,izb) +
     &         (tTwo*Rnxyz(iZeta,1,ixa+1,ixb)*Alpha(iZeta) +
     &            xa*Rnxyz(iZeta,1,ixa-1,ixb) )*
     &               Rnxyz(iZeta,2,iya,iyb)*
     &               Txyz(iZeta,3,iza,izb) )
 30               Continue
               Else
                  Do 31 iZeta = 1, nZeta
                     Final(iZeta,ipa,ipb,1) =
     &                rKappa(iZeta) * Alpha(iZeta) * (
     &           tTwo*Txyz(iZeta,1,ixa+1,ixb) *
     &               Rnxyz(iZeta,2,iya,iyb)*
     &               Rnxyz(iZeta,3,iza,izb) +
     &          tTwo*Rnxyz(iZeta,1,ixa+1,ixb) *
     &               Txyz(iZeta,2,iya,iyb)*
     &               Rnxyz(iZeta,3,iza,izb) +
     &          tTwo*Rnxyz(iZeta,1,ixa+1,ixb) *
     &               Rnxyz(iZeta,2,iya,iyb)*
     &               Txyz(iZeta,3,iza,izb) )
 31               Continue
               End If
            End If
            If (IfGrad(1,2)) Then
               If (ixb.gt.0) Then
                  xb = Dble(-ixb)
                  Do 35 iZeta = 1, nZeta
                     Final(iZeta,ipa,ipb,1) =
     &                rKappa(iZeta) * (
     &          (tTwo*Txyz(iZeta,1,ixa,ixb+1)*Beta(iZeta) +
     &             xb*Txyz(iZeta,1,ixa,ixb-1) )*
     &               Rnxyz(iZeta,2,iya,iyb)*
     &               Rnxyz(iZeta,3,iza,izb) +
     &         (tTwo*Rnxyz(iZeta,1,ixa,ixb+1)*Beta(iZeta) +
     &            xb*Rnxyz(iZeta,1,ixa,ixb-1) )*
     &               Txyz(iZeta,2,iya,iyb)*
     &               Rnxyz(iZeta,3,iza,izb) +
     &         (tTwo*Rnxyz(iZeta,1,ixa,ixb+1)*Beta(iZeta) +
     &            xb*Rnxyz(iZeta,1,ixa,ixb-1) )*
     &               Rnxyz(iZeta,2,iya,iyb)*
     &               Txyz(iZeta,3,iza,izb) )
 35               Continue
               Else
                  Do 36 iZeta = 1, nZeta
                     Final(iZeta,ipa,ipb,1) =
     &                rKappa(iZeta) * Beta(iZeta) * (
     &           tTwo*Txyz(iZeta,1,ixa,ixb+1) *
     &               Rnxyz(iZeta,2,iya,iyb)*
     &               Rnxyz(iZeta,3,iza,izb) +
     &          tTwo*Rnxyz(iZeta,1,ixa,ixb+1) *
     &               Txyz(iZeta,2,iya,iyb)*
     &               Rnxyz(iZeta,3,iza,izb) +
     &          tTwo*Rnxyz(iZeta,1,ixa,ixb+1) *
     &               Rnxyz(iZeta,2,iya,iyb)*
     &               Txyz(iZeta,3,iza,izb) )
 36               Continue
               End If
            End If
            If (IfGrad(2,1)) Then
               If (iya.gt.0) Then
                  ya = Dble(-iya)
                  Do 40 iZeta = 1, nZeta
                     Final(iZeta,ipa,ipb,1) =
     &                rKappa(iZeta) * (
     &                Txyz(iZeta,1,ixa,ixb)*
     &         (tTwo*Rnxyz(iZeta,2,iya+1,iyb)*Alpha(iZeta) +
     &            ya*Rnxyz(iZeta,2,iya-1,iyb) )*
     &               Rnxyz(iZeta,3,iza,izb) +
     &               Rnxyz(iZeta,1,ixa,ixb)*
     &          (tTwo*Txyz(iZeta,2,iya+1,iyb)*Alpha(iZeta) +
     &             ya*Txyz(iZeta,2,iya-1,iyb) )*
     &               Rnxyz(iZeta,3,iza,izb) +
     &               Rnxyz(iZeta,1,ixa,ixb) *
     &         (tTwo*Rnxyz(iZeta,2,iya+1,iyb)*Alpha(iZeta) +
     &            ya*Rnxyz(iZeta,2,iya-1,iyb) )*
     &               Txyz(iZeta,3,iza,izb) )
 40               Continue
               Else
                  Do 41 iZeta = 1, nZeta
                     Final(iZeta,ipa,ipb,1) =
     &                rKappa(iZeta) * Alpha(iZeta) * (
     &                Txyz(iZeta,1,ixa,ixb)*
     &          tTwo*Rnxyz(iZeta,2,iya+1,iyb) *
     &               Rnxyz(iZeta,3,iza,izb) +
     &               Rnxyz(iZeta,1,ixa,ixb)*
     &           tTwo*Txyz(iZeta,2,iya+1,iyb) *
     &               Rnxyz(iZeta,3,iza,izb) +
     &               Rnxyz(iZeta,1,ixa,ixb) *
     &          tTwo*Rnxyz(iZeta,2,iya+1,iyb) *
     &               Txyz(iZeta,3,iza,izb) )
 41               Continue
               End If
            End If
            If (IfGrad(2,2)) Then
               If (iyb.gt.0) Then
                  yb = Dble(-iyb)
                  Do 45 iZeta = 1, nZeta
                     Final(iZeta,ipa,ipb,1) =
     &                rKappa(iZeta) * (
     &                Txyz(iZeta,1,ixa,ixb)*
     &         (tTwo*Rnxyz(iZeta,2,iya,iyb+1)*Beta(iZeta) +
     &            yb*Rnxyz(iZeta,2,iya,iyb-1) )*
     &               Rnxyz(iZeta,3,iza,izb) +
     &               Rnxyz(iZeta,1,ixa,ixb)*
     &          (tTwo*Txyz(iZeta,2,iya,iyb+1)*Beta(iZeta) +
     &             yb*Txyz(iZeta,2,iya,iyb-1) )*
     &               Rnxyz(iZeta,3,iza,izb) +
     &               Rnxyz(iZeta,1,ixa,ixb) *
     &         (tTwo*Rnxyz(iZeta,2,iya,iyb+1)*Beta(iZeta) +
     &            yb*Rnxyz(iZeta,2,iya,iyb-1) )*
     &               Txyz(iZeta,3,iza,izb) )
 45               Continue
               Else
                  Do 46 iZeta = 1, nZeta
                     Final(iZeta,ipa,ipb,1) =
     &                rKappa(iZeta) * Beta(iZeta) * (
     &                Txyz(iZeta,1,ixa,ixb)*
     &          tTwo*Rnxyz(iZeta,2,iya,iyb+1) *
     &               Rnxyz(iZeta,3,iza,izb) +
     &               Rnxyz(iZeta,1,ixa,ixb)*
     &           tTwo*Txyz(iZeta,2,iya,iyb+1) *
     &               Rnxyz(iZeta,3,iza,izb) +
     &               Rnxyz(iZeta,1,ixa,ixb) *
     &          tTwo*Rnxyz(iZeta,2,iya,iyb+1) *
     &               Txyz(iZeta,3,iza,izb) )
 46               Continue
               End If
            End If
            If (IfGrad(3,1)) Then
               If (iza.gt.0) Then
                  za = Dble(-iza)
                  Do 50 iZeta = 1, nZeta
                     Final(iZeta,ipa,ipb,1) =
     &                rKappa(iZeta) * (
     &                Txyz(iZeta,1,ixa,ixb)*
     &               Rnxyz(iZeta,2,iya,iyb)*
     &         (tTwo*Rnxyz(iZeta,3,iza+1,izb)*Alpha(iZeta) +
     &            za*Rnxyz(iZeta,3,iza-1,izb) ) +
     &               Rnxyz(iZeta,1,ixa,ixb)*
     &                Txyz(iZeta,2,iya,iyb)*
     &         (tTwo*Rnxyz(iZeta,3,iza+1,izb)*Alpha(iZeta) +
     &            za*Rnxyz(iZeta,3,iza-1,izb) ) +
     &               Rnxyz(iZeta,1,ixa,ixb) *
     &               Rnxyz(iZeta,2,iya,iyb)*
     &          (tTwo*Txyz(iZeta,3,iza+1,izb)*Alpha(iZeta) +
     &             za*Txyz(iZeta,3,iza-1,izb) ) )
 50               Continue
               Else
                  Do 51 iZeta = 1, nZeta
                     Final(iZeta,ipa,ipb,1) =
     &                rKappa(iZeta) * Alpha(iZeta) * (
     &                Txyz(iZeta,1,ixa,ixb)*
     &               Rnxyz(iZeta,2,iya,iyb)*
     &          tTwo*Rnxyz(iZeta,3,iza+1,izb) +
     &               Rnxyz(iZeta,1,ixa,ixb)*
     &                Txyz(iZeta,2,iya,iyb)*
     &          tTwo*Rnxyz(iZeta,3,iza+1,izb) +
     &               Rnxyz(iZeta,1,ixa,ixb) *
     &               Rnxyz(iZeta,2,iya,iyb)*
     &           tTwo*Txyz(iZeta,3,iza+1,izb) )
 51               Continue
               End If
            End If
            If (IfGrad(3,2)) Then
               If (izb.gt.0) Then
                  zb = Dble(-izb)
                  Do 55 iZeta = 1, nZeta
                     Final(iZeta,ipa,ipb,1) =
     &                rKappa(iZeta) * (
     &                Txyz(iZeta,1,ixa,ixb)*
     &               Rnxyz(iZeta,2,iya,iyb)*
     &         (tTwo*Rnxyz(iZeta,3,iza,izb+1)*Beta(iZeta) +
     &            zb*Rnxyz(iZeta,3,iza,izb-1) ) +
     &               Rnxyz(iZeta,1,ixa,ixb)*
     &                Txyz(iZeta,2,iya,iyb)*
     &         (tTwo*Rnxyz(iZeta,3,iza,izb+1)*Beta(iZeta) +
     &            zb*Rnxyz(iZeta,3,iza,izb-1) ) +
     &               Rnxyz(iZeta,1,ixa,ixb) *
     &               Rnxyz(iZeta,2,iya,iyb)*
     &          (tTwo*Txyz(iZeta,3,iza,izb+1)*Beta(iZeta) +
     &             zb*Txyz(iZeta,3,iza,izb-1) ) )
 55               Continue
               Else
                  Do 56 iZeta = 1, nZeta
                     Final(iZeta,ipa,ipb,1) =
     &                rKappa(iZeta) * Beta(iZeta) * (
     &                Txyz(iZeta,1,ixa,ixb)*
     &               Rnxyz(iZeta,2,iya,iyb)*
     &          tTwo*Rnxyz(iZeta,3,iza,izb+1) +
     &               Rnxyz(iZeta,1,ixa,ixb)*
     &                Txyz(iZeta,2,iya,iyb)*
     &          tTwo*Rnxyz(iZeta,3,iza,izb+1) +
     &               Rnxyz(iZeta,1,ixa,ixb) *
     &               Rnxyz(iZeta,2,iya,iyb)*
     &           tTwo*Txyz(iZeta,3,iza,izb+1) )
 56               Continue
               End If
            End If
*
 21      Continue
 20      Continue
 11   Continue
 10   Continue
*
*
c     Call QExit('CmbnT1')
      Return
      End
