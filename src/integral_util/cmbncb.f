!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!                                                                      *
! Copyright (C) 1991, Roland Lindh                                     *
!***********************************************************************
!#define _DEBUGPRINT_
      SubRoutine CmbnCB(Rnxyz,nZeta,la,lb,rKappa,Final,Beta,
     &                  IfGrad,ld,nVecCB)
!***********************************************************************
!                                                                      *
! Object: compute the gradient of the overlap matrix.                  *
!                                                                      *
!     Author: Roland Lindh, Dept. of Theoretical Chemistry,            *
!             University of Lund, SWEDEN                               *
!             October '91.                                             *
!***********************************************************************
      use Constants, only: Two
      Implicit None
      Integer nZeta, la, lb, ld, nVecCB
      Real*8 Final(nZeta,(la+1)*(la+2)/2,(lb+1)*(lb+2)/2,4),
     &       rKappa(nZeta), Beta(nZeta),
     &       Rnxyz(nZeta,3,0:la,0:lb+ld)
      Logical IfGrad(3)

      Integer ixyz, ix, iz, Ind
      Integer ixa, ixb, iya, iyb, iza, izb, ipa, ipb, iZeta,
     &        iyaMax, iybMax
      Real*8 xB, yB, zB, tTwo
!
!     Statement function for Cartesian index
!
      Ind(ixyz,ix,iz) = (ixyz-ix)*(ixyz-ix+1)/2 + iz + 1
!
#ifdef _DEBUGPRINT_
      Call RecPrt(' In CmbnCB: rKappa',' ',rKappa,1,nZeta)
      Call RecPrt(' In CmbnCB: Beta  ',' ',Beta  ,1,nZeta)
#endif
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
!
!           Combine overlap integral gradients
!
            nVecCB = 1
            Do 26 iZeta = 1, nZeta
               Final(iZeta,ipa,ipb,nVecCB) = rKappa(iZeta)*
     &                 Rnxyz(iZeta,1,ixa,ixb)*
     &                 Rnxyz(iZeta,2,iya,iyb)*
     &                 Rnxyz(iZeta,3,iza,izb)
 26         Continue
            tTwo = Two
            If (IfGrad(1)) Then
               nVecCB = nVecCB+1
               If (ixb.gt.0) Then
                  xb = DBLE(-ixb)
                  Do 35 iZeta = 1, nZeta
                     Final(iZeta,ipa,ipb,nVecCB) = rKappa(iZeta)*
     &     (tTwo*Beta(iZeta)*Rnxyz(iZeta,1,ixa,ixb+1) +
     &                    xb*Rnxyz(iZeta,1,ixa,ixb-1))*
     &                       Rnxyz(iZeta,2,iya,iyb)*
     &                       Rnxyz(iZeta,3,iza,izb)
 35               Continue
               Else
                  Do 36 iZeta = 1, nZeta
                     Final(iZeta,ipa,ipb,nVecCB) = rKappa(iZeta)*
     &      tTwo*Beta(iZeta)*Rnxyz(iZeta,1,ixa,ixb+1)*
     &                       Rnxyz(iZeta,2,iya,iyb)*
     &                       Rnxyz(iZeta,3,iza,izb)
 36               Continue
               End If
            End If
            If (IfGrad(2)) Then
               nVecCB = nVecCB+1
               If (iyb.gt.0) Then
                  yb = DBLE(-iyb)
                  Do 45 iZeta = 1, nZeta
                     Final(iZeta,ipa,ipb,nVecCB) = rKappa(iZeta)*
     &                       Rnxyz(iZeta,1,ixa,ixb)*
     &     (tTwo*Beta(iZeta)*Rnxyz(iZeta,2,iya,iyb+1) +
     &                    yb*Rnxyz(iZeta,2,iya,iyb-1))*
     &                       Rnxyz(iZeta,3,iza,izb)
 45               Continue
               Else
                  Do 46 iZeta = 1, nZeta
                     Final(iZeta,ipa,ipb,nVecCB) = rKappa(iZeta)*
     &                       Rnxyz(iZeta,1,ixa,ixb)*
     &      tTwo*Beta(iZeta)*Rnxyz(iZeta,2,iya,iyb+1)*
     &                       Rnxyz(iZeta,3,iza,izb)
 46               Continue
               End If
            End If
            If (IfGrad(3)) Then
               nVecCB = nVecCB+1
               If (izb.gt.0) Then
                  zb = DBLE(-izb)
                  Do 55 iZeta = 1, nZeta
                     Final(iZeta,ipa,ipb,nVecCB) = rKappa(iZeta)*
     &                       Rnxyz(iZeta,1,ixa,ixb)*
     &                       Rnxyz(iZeta,2,iya,iyb)*
     &     (tTwo*Beta(iZeta)*Rnxyz(iZeta,3,iza,izb+1) +
     &                    zb*Rnxyz(iZeta,3,iza,izb-1))
 55               Continue
               Else
                  Do 56 iZeta = 1, nZeta
                     Final(iZeta,ipa,ipb,nVecCB) = rKappa(iZeta)*
     &                       Rnxyz(iZeta,1,ixa,ixb)*
     &                       Rnxyz(iZeta,2,iya,iyb)*
     &      tTwo*Beta(iZeta)*Rnxyz(iZeta,3,iza,izb+1)
 56               Continue
               End If
            End If
!
 21      Continue
 20      Continue
 11   Continue
 10   Continue
!
      Return
      End SubRoutine CmbnCB
