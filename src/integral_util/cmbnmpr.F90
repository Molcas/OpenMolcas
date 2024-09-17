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
subroutine CmbnMPr(Rnr,nZeta,la,lb,lr,rFinal,nComp)
!***********************************************************************
!     Author: K.Pfingst                                                *
!***********************************************************************

use rmat, only: GammaPh, GammaTh
use Definitions, only: wp, iwp
#ifdef _DEBUGPRINT_
use Definitions, only: u6
#endif

implicit none
integer(kind=iwp), intent(in) :: nZeta, la, lb, lr, nComp
real(kind=wp), intent(in) :: Rnr(nZeta,0:la+lb+lr)
real(kind=wp), intent(out) :: rFinal(nZeta,nComp,(la+1)*(la+2)/2,(lb+1)*(lb+2)/2)
integer(kind=iwp) :: iComp, ipa, ipb, ixa, ixb, iya, iyaMax, iyb, iybMax, iz, iza, izb, iZeta, lcosf, lcost, lrs, lsinf, lsint
real(kind=wp) :: Fact
! Statement function for Cartesian index
integer(kind=iwp) :: ixyz, ix, iy, Ind
Ind(ixyz,ix,iz) = (ixyz-ix)*(ixyz-ix+1)/2+iz+1

do ixa=0,la
  iyaMax = la-ixa
  do ixb=0,lb
    iybMax = lb-ixb
    do iya=0,iyaMax
      iza = la-ixa-iya
      ipa = Ind(la,ixa,iza)
      do iyb=0,iybMax
        izb = lb-ixb-iyb
        ipb = Ind(lb,ixb,izb)
#       ifdef _DEBUGPRINT_
        write(u6,*) ixa,iya,iza,ixb,iyb,izb
        write(u6,*) ipa,ipb
#       endif

        ! Combine multipole moment integrals

        iComp = 0
        do ix=lr,0,-1
          do iy=lr-ix,0,-1
            iz = lr-ix-iy
            iComp = iComp+1
            lrs = ixa+ixb+ix+iya+iyb+iy+iza+izb+iz
            lcost = iza+izb+iz
            lsint = ixa+ixb+ix+iya+iyb+iy
            lsinf = iya+iyb+iy
            lcosf = ixa+ixb+ix
            Fact = gammath(lsint,lcost)*gammaph(lsinf,lcosf)
            do iZeta=1,nZeta
              rFinal(iZeta,iComp,ipa,ipb) = Fact*Rnr(iZeta,lrs)
            end do
          end do
        end do

      end do
    end do
  end do
end do

return

end subroutine CmbnMPr
