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
! Copyright (C) 1991,1992, Roland Lindh                                *
!***********************************************************************

!#define _DEBUGPRINT_
subroutine CmbnRF(Rnxyz,nZeta,la,lb,lr,Zeta,rKappa,rFinal,nComp,Fact,Temp)
!***********************************************************************
!     Author: Roland Lindh, Dept. of Theoretical Chemistry,            *
!             University of Lund, SWEDEN                               *
!             Modified for reaction field calculations July '92        *
!***********************************************************************

use Index_Functions, only: C_Ind, nTri_Elem1, nTri3_Elem
use Definitions, only: wp, iwp
#ifdef _DEBUGPRINT_
use Definitions, only: u6
#endif

implicit none
integer(kind=iwp), intent(in) :: nZeta, la, lb, lr, nComp
real(kind=wp), intent(in) :: Rnxyz(nZeta,3,0:la,0:lb,0:lr), Zeta(nZeta), rKappa(nZeta)
real(kind=wp), intent(out) :: rFinal(nZeta,nTri_Elem1(la),nTri_Elem1(lb),nComp), Fact(nZeta), Temp(nZeta)
integer(kind=iwp) :: iComp, ipa, ipb, ir, ix, ixa, ixb, iy, iya, iyaMax, iyb, iybMax, iz, iza, izb

Fact(:) = rKappa(:)*sqrt(Zeta(:)**(-3))
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

        do ix=0,lr
          do iy=0,lr-ix
            Temp(:) = Fact(:)*Rnxyz(:,1,ixa,ixb,ix)*Rnxyz(:,2,iya,iyb,iy)
            do ir=ix+iy,lr
              iz = ir-ix-iy
              iComp = C_Ind(ir,ix,iz)+nTri3_Elem(ir)
              rFinal(:,ipa,ipb,iComp) = Temp(:)*Rnxyz(:,3,iza,izb,iz)
            end do
          end do
        end do

      end do
    end do
  end do
end do

#ifdef _DEBUGPRINT_
call RecPrt('Final',' ',rFinal,nZeta*nTri_Elem1(la)*nTri_Elem1(lb),nComp)
#endif

return

end subroutine CmbnRF
