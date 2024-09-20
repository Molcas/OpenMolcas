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

use Index_Functions, only: C_Ind, nTri_Elem1
use rmat, only: GammaPh, GammaTh
use Definitions, only: wp, iwp
#ifdef _DEBUGPRINT_
use Definitions, only: u6
#endif

implicit none
integer(kind=iwp), intent(in) :: nZeta, la, lb, lr, nComp
real(kind=wp), intent(in) :: Rnr(nZeta,0:la+lb+lr)
real(kind=wp), intent(out) :: rFinal(nZeta,nComp,nTri_Elem1(la),nTri_Elem1(lb))
integer(kind=iwp) :: iComp, ipa, ipb, ix, ixa, ixb, iy, iya, iyaMax, iyb, iybMax, iz, iza, izb, lcosf, lcost, lrs, lsinf, lsint
real(kind=wp) :: Fact

do ixa=0,la
  iyaMax = la-ixa
  do ixb=0,lb
    iybMax = lb-ixb
    do iya=0,iyaMax
      iza = la-ixa-iya
      ipa = C_Ind(la,ixa,iza)
      do iyb=0,iybMax
        izb = lb-ixb-iyb
        ipb = C_Ind(lb,ixb,izb)
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
            rFinal(:,iComp,ipa,ipb) = Fact*Rnr(:,lrs)
          end do
        end do

      end do
    end do
  end do
end do

return

end subroutine CmbnMPr
