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

subroutine CmbnMP(Rnxyz,nZeta,la,lb,lr,Zeta,rKappa,rFinal,nComp)
!***********************************************************************
!     Author: Roland Lindh, Dept. of Theoretical Chemistry,            *
!             University of Lund, SWEDEN                               *
!***********************************************************************

use Index_Functions, only: C_Ind, nTri_Elem1
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: nZeta, la, lb, lr, nComp
real(kind=wp), intent(in) :: Rnxyz(nZeta,3,0:la,0:lb,0:lr), Zeta(nZeta), rKappa(nZeta)
real(kind=wp), intent(inout) :: rFinal(nZeta,nTri_Elem1(la),nTri_Elem1(lb),nComp)
integer(kind=iwp) :: iComp, ipa, ipb, ix, ixa, ixb, iy, iya, iyaMax, iyb, iybMax, iz, iza, izb

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
        !if (iPrint >= 99) then
        !  write(u6,*) ixa,iya,iza,ixb,iyb,izb
        !  write(u6,*) ipa,ipb
        !end if

        ! Combine multipole moment integrals

        iComp = 0
        do ix=lr,0,-1
          do iy=lr-ix,0,-1
            iz = lr-ix-iy
            iComp = iComp+1
            !write(u6,*) ix,iy,iz,iComp
            rFinal(:,ipa,ipb,iComp) = rKappa(:)/sqrt(Zeta(:)**3)*Rnxyz(:,1,ixa,ixb,ix)*Rnxyz(:,2,iya,iyb,iy)*Rnxyz(:,3,iza,izb,iz)
          end do
        end do

      end do
    end do
  end do
end do

return

end subroutine CmbnMP
