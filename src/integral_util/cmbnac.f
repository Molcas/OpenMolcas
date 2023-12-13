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
      SubRoutine CmbnAC(Rnxyz,nZeta,la,lb,rKappa,Final,Alpha,
     &                  IfGrad,ld,nVecAC)
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
      Integer nZeta, la, lb, ld, nVecAC
      Real*8 Final(nZeta,(la+1)*(la+2)/2,(lb+1)*(lb+2)/2,4),
     &       rKappa(nZeta),
     &       Rnxyz(nZeta,3,0:la+ld,0:lb), Alpha(nZeta)
      Logical IfGrad(3)

      Integer ixa, ixb, iya, iyb, iza, izb, iZeta, iyaMax, iybMax,
     &        ipa, ipb
      Integer ixyz, ix, iz, Ind
      Real*8 tTwo, XA, YA, ZA
!
!     Statement function for Cartesian index
!
      Ind(ixyz,ix,iz) = (ixyz-ix)*(ixyz-ix+1)/2 + iz + 1
!
!
#ifdef _DEBUGPRINT_
      Call RecPrt(' In CmbnAC: rKappa',' ',rKappa,1,nZeta)
      Call RecPrt(' In CmbnAC: Alpha ',' ',Alpha ,1,nZeta)
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
            nVecAC = 1
            Do 35 iZeta = 1, nZeta
               Final(iZeta,ipa,ipb,nVecAC) = rKappa(iZeta)*
     &                 Rnxyz(iZeta,1,ixa,ixb)*
     &                 Rnxyz(iZeta,2,iya,iyb)*
     &                 Rnxyz(iZeta,3,iza,izb)
 35         Continue
            tTwo = Two
            If (IfGrad(1)) Then
               nVecAC = nVecAC + 1
               If (ixa.gt.0) Then
                  xa = DBLE(-ixa)
                  Do 30 iZeta = 1, nZeta
                     Final(iZeta,ipa,ipb,nVecAC) = rKappa(iZeta)*
     &    (tTwo*Alpha(iZeta)*Rnxyz(iZeta,1,ixa+1,ixb) +
     &                    xa*Rnxyz(iZeta,1,ixa-1,ixb))*
     &                       Rnxyz(iZeta,2,iya,iyb)*
     &                       Rnxyz(iZeta,3,iza,izb)
 30               Continue
               Else
                  Do 31 iZeta = 1, nZeta
                     Final(iZeta,ipa,ipb,nVecAC) = rKappa(iZeta)*
     &     tTwo*Alpha(iZeta)*Rnxyz(iZeta,1,ixa+1,ixb)*
     &                       Rnxyz(iZeta,2,iya,iyb)*
     &                       Rnxyz(iZeta,3,iza,izb)
 31               Continue
               End If
            End If
            If (IfGrad(2)) Then
               nVecAC = nVecAC + 1
               If (iya.gt.0) Then
                  ya = DBLE(-iya)
                  Do 40 iZeta = 1, nZeta
                     Final(iZeta,ipa,ipb,nVecAC) = rKappa(iZeta)*
     &                       Rnxyz(iZeta,1,ixa,ixb)*
     &    (tTwo*Alpha(iZeta)*Rnxyz(iZeta,2,iya+1,iyb) +
     &                    ya*Rnxyz(iZeta,2,iya-1,iyb))*
     &                       Rnxyz(iZeta,3,iza,izb)
 40               Continue
               Else
                  Do 41 iZeta = 1, nZeta
                     Final(iZeta,ipa,ipb,nVecAC) = rKappa(iZeta)*
     &                       Rnxyz(iZeta,1,ixa,ixb)*
     &     tTwo*Alpha(iZeta)*Rnxyz(iZeta,2,iya+1,iyb)*
     &                       Rnxyz(iZeta,3,iza,izb)
 41               Continue
               End If
            End If
            If (IfGrad(3)) Then
               nVecAC = nVecAC + 1
               If (iza.gt.0) Then
                  za = DBLE(-iza)
                  Do 50 iZeta = 1, nZeta
                     Final(iZeta,ipa,ipb,nVecAC) = rKappa(iZeta)*
     &                       Rnxyz(iZeta,1,ixa,ixb)*
     &                       Rnxyz(iZeta,2,iya,iyb)*
     &    (tTwo*Alpha(iZeta)*Rnxyz(iZeta,3,iza+1,izb) +
     &                    za*Rnxyz(iZeta,3,iza-1,izb))
 50               Continue
               Else
                  Do 51 iZeta = 1, nZeta
                     Final(iZeta,ipa,ipb,nVecAC) = rKappa(iZeta)*
     &                       Rnxyz(iZeta,1,ixa,ixb)*
     &                       Rnxyz(iZeta,2,iya,iyb)*
     &     tTwo*Alpha(iZeta)*Rnxyz(iZeta,3,iza+1,izb)
 51               Continue
               End If
            End If
!
 21      Continue
 20      Continue
 11   Continue
 10   Continue
!
      Return
      End SubRoutine CmbnAC
