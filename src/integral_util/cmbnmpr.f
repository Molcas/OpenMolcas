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
! Copyright (C) Kurt Pfingst                                           *
!***********************************************************************
!#define _DEBUGPRINT_
      SubRoutine CmbnMPr(Rnr,nZeta,la,lb,lr,Final,nComp)
!***********************************************************************
!     Author: K.Pfingst                                                *
!***********************************************************************
      use rmat, only: lCosT, lSinT, lSinF, lCosF, GammaPh, GammaTh
      Implicit None
      Integer nZeta, nComp, la, lb, lr
      Real*8 Final(nZeta,nComp,(la+1)*(la+2)/2,(lb+1)*(lb+2)/2),
     &       Rnr(nZeta,0:(la+lb+lr))
      Integer ixa, ixb, iya, iyb, iza, izb, iyaMax, iybMax, ipa, ipb,
     &        iComp, lrs, iZeta, iz
      Real*8 Fact

      Integer ixyz, ix, iy, Ind
!
!     Statement function for Cartesian index
!
      Ind(ixyz,ix,iz) = (ixyz-ix)*(ixyz-ix+1)/2 + iz + 1
!
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
#ifdef _DEBUGPRINT_
            Write (6,*) ixa,iya,iza,ixb,iyb,izb
            Write (6,*) ipa,ipb
#endif
!
!           Combine multipole moment integrals
!
            iComp = 0
            Do 41 ix = lr, 0, -1
               Do 42 iy = lr-ix, 0, -1
                  iz = lr-ix-iy
                  iComp=iComp+1
                  lrs=ixa+ixb+ix+iya+iyb+iy+iza+izb+iz
                  lcost=iza+izb+iz
                  lsint=ixa+ixb+ix+iya+iyb+iy
                  lsinf=iya+iyb+iy
                  lcosf=ixa+ixb+ix
                  Fact=gammath(lsint,lcost)*gammaph(lsinf,lcosf)
                  Do 30 iZeta = 1, nZeta
                     Final(iZeta,iComp,ipa,ipb) = Fact *
     *                       Rnr(iZeta,lrs)
30                Continue
42             Continue
41          Continue
!
21       Continue
20       Continue
11    Continue
10    Continue
!
      Return
      End SubRoutine CmbnMPr
