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

subroutine CCmbnMP(Rnxyz,nZeta,la,lb,lr,Zeta,rKappa,rFinal,nComp,kVector,P)
!***********************************************************************
!     Author: Roland Lindh, Dept. of Theoretical Chemistry,            *
!             University of Lund, SWEDEN                               *
!***********************************************************************

use Index_Functions, only: C_Ind, nTri_Elem1
use Constants, only: One, Quart, Onei
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: nZeta, la, lb, lr, nComp
complex(kind=wp), intent(in) :: Rnxyz(nZeta,3,0:la,0:lb,0:lr)
real(kind=wp), intent(in) :: Zeta(nZeta), rKappa(nZeta), kVector(3), P(nZeta,3)
real(kind=wp), intent(out) :: rFinal(nZeta,nTri_Elem1(la),nTri_Elem1(lb),nComp)
integer(kind=iwp) :: iComp, ipa, ipb, ix, ixa, ixb, iy, iya, iyaMax, iyb, iybMax, iz, iza, izb, iZeta
complex(kind=wp) :: Temp
real(kind=wp) :: Fact, k_Dot_P, kModQ, rTemp

kModQ = Quart*(kVector(1)**2+kVector(2)**2+kVector(3)**2)

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

        ! Combine multipole moment integrals

        iComp = 0
        do ix=lr,0,-1
          do iy=lr-ix,0,-1
            iz = lr-ix-iy
            do iZeta=1,nZeta
              rTemp = kModQ/Zeta(iZeta)
              Fact = rKappa(iZeta)*(One/sqrt(Zeta(iZeta)**3))*exp(-rTemp)
              k_Dot_P = kVector(1)*P(iZeta,1)+kVector(2)*P(iZeta,2)+kVector(3)*P(iZeta,3)
              Temp = exp(Onei*k_Dot_P)*Fact*Rnxyz(iZeta,1,ixa,ixb,ix)*Rnxyz(iZeta,2,iya,iyb,iy)*Rnxyz(iZeta,3,iza,izb,iz)
              rFinal(iZeta,ipa,ipb,iComp+1) = real(Temp)
              rFinal(iZeta,ipa,ipb,iComp+2) = aimag(Temp)
            end do
            iComp = iComp+2
          end do
        end do

      end do
    end do
  end do
end do

return

end subroutine CCmbnMP
