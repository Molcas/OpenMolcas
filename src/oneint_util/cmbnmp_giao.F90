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
! Copyright (C) 1991,2000, Roland Lindh                                *
!***********************************************************************

subroutine CmbnMP_GIAO(Rnxyz,nZeta,la,lb,lr,Zeta,rKappa,rFinal,nComp,nB,RAB,C)
!***********************************************************************
!     Author: Roland Lindh, Dept. of Theoretical Chemistry,            *
!             University of Lund, SWEDEN                               *
!                                                                      *
!             Modified to GIAO 1st derivatives by R. Lindh in Tokyo,   *
!              Japan, January 2000.                                    *
!***********************************************************************

use Index_Functions, only: C_Ind, nTri_Elem1
use Constants, only: Half, OneHalf
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: nZeta, la, lb, lr, nComp, nB
real(kind=wp), intent(in) :: Rnxyz(nZeta,3,0:la,0:lb,0:lr+1), Zeta(nZeta), rKappa(nZeta), RAB(3), C(3)
real(kind=wp), intent(out) :: rFinal(nZeta,nTri_Elem1(la),nTri_Elem1(lb),nComp,nB)
integer(kind=iwp) :: iBx, iBy, iBz, iComp, ipa, ipb, ix, ix_(3,2), ixa, ixb, iy, iya, iyaMax, iyb, iybMax, iz, iza, izb, iZeta
real(kind=wp) :: Fact, temp, tempy, tempz

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

        do iBx=1,3
          iBy = iBx+1
          if (iBy > 3) iBy = iBy-3
          iBz = iBy+1
          if (iBz > 3) iBz = iBz-3
          ix_(:,:) = 0
          ix_(iBz,:) = 1
          iComp = 0
          do ix=lr,0,-1
            do iy=lr-ix,0,-1
              iz = lr-ix-iy
              iComp = iComp+1
              do iZeta=1,nZeta
                Fact = rKappa(iZeta)*Zeta(iZeta)**(-OneHalf)
                temp = Rnxyz(iZeta,1,ixa,ixb,ix)*Rnxyz(iZeta,2,iya,iyb,iy)*Rnxyz(iZeta,3,iza,izb,iz)
                tempz = Rnxyz(iZeta,1,ixa,ixb,ix+ix_(1,1))*Rnxyz(iZeta,2,iya,iyb,iy+ix_(2,1))*Rnxyz(iZeta,3,iza,izb,iz+ix_(3,1))
                tempy = Rnxyz(iZeta,1,ixa,ixb,ix+ix_(1,2))*Rnxyz(iZeta,2,iya,iyb,iy+ix_(2,2))*Rnxyz(iZeta,3,iza,izb,iz+ix_(3,2))

                ! The term has only an imaginary component

                rFinal(iZeta,ipa,ipb,iComp,iBx) = Half*Fact*(RAB(iBy)*(Tempz+C(iBz)*Temp)-RAB(iBz)*(Tempy+C(iBy)*Temp))
              end do
            end do
          end do
        end do  ! iB

      end do
    end do
  end do
end do

return

end subroutine CmbnMP_GIAO
