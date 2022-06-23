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

subroutine CmbnMV(Rnxyz,nZeta,la,lb,lr,Zeta,rKappa,rFinal,nComp,rV2Int,rV4Int)
!***********************************************************************
!     Author: Roland Lindh, Dept. of Theoretical Chemistry,            *
!             University of Lund, SWEDEN, February '91                 *
!***********************************************************************

use Index_Functions, only: nTri_Elem1
use Constants, only: One, Two, Three, Eight, c_in_au
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp) :: nZeta, la, lb, lr, nComp
real(kind=wp) :: Rnxyz(nZeta,3,0:la+2,0:lb+2,0:lr), Zeta(nZeta), rKappa(nZeta), rFinal(nZeta,nComp,nTri_Elem1(la),nTri_Elem1(lb)), &
                 rV2Int(nZeta,3,0:la,0:lb,2), rV4Int(nZeta,3,0:la,0:lb)
integer(kind=iwp) :: iComp, ipa, ipb, ixa, ixb, iya, iyaMax, iyb, iybMax, iza, izb, iZeta
real(kind=wp) :: Fact, rMVel, x2x2, x2y2, x2z2, y2x2, y2y2, y2z2, z2x2, z2y2, z2z2
real(kind=wp), parameter :: Const = -One/(Eight*c_in_au**2)
! Statement function for Cartesian index
integer(kind=iwp) :: Ind, ixyz, ix, iz
Ind(ixyz,ix,iz) = (ixyz-ix)*(ixyz-ix+1)/2+iz+1

!iRout = 191
!iPrint = nPrint(iRout)

iComp = 1
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
        !if (iPrint >= 99) then
        !  write(u6,*) ixa,iya,iza,ixb,iyb,izb
        !  write(u6,*) ipa,ipb
        !end if

        ! Combine integrals

        do iZeta=1,nZeta
          Fact = rKappa(iZeta)*Zeta(iZeta)**(-Three/Two)*Const
          x2x2 = rV4Int(iZeta,1,ixa,ixb)*Rnxyz(iZeta,2,iya,iyb,0)*Rnxyz(iZeta,3,iza,izb,0)
          x2y2 = rV2Int(iZeta,1,ixa,ixb,1)*rV2Int(iZeta,2,iya,iyb,2)*Rnxyz(iZeta,3,iza,izb,0)
          x2z2 = rV2Int(iZeta,1,ixa,ixb,1)*Rnxyz(iZeta,2,iya,iyb,0)*rV2Int(iZeta,3,iza,izb,2)
          y2x2 = rV2Int(iZeta,1,ixa,ixb,2)*rV2Int(iZeta,2,iya,iyb,1)*Rnxyz(iZeta,3,iza,izb,0)
          y2y2 = Rnxyz(iZeta,1,ixa,ixb,0)*rV4Int(iZeta,2,iya,iyb)*Rnxyz(iZeta,3,iza,izb,0)
          y2z2 = Rnxyz(iZeta,1,ixa,ixb,0)*rV2Int(iZeta,2,iya,iyb,1)*rV2Int(iZeta,3,iza,izb,2)
          z2x2 = rV2Int(iZeta,1,ixa,ixb,2)*Rnxyz(iZeta,2,iya,iyb,0)*rV2Int(iZeta,3,iza,izb,1)
          z2y2 = Rnxyz(iZeta,1,ixa,ixb,0)*rV2Int(iZeta,2,iya,iyb,2)*rV2Int(iZeta,3,iza,izb,1)
          z2z2 = Rnxyz(iZeta,1,ixa,ixb,0)*Rnxyz(iZeta,2,iya,iyb,0)*rV4Int(iZeta,3,iza,izb)

          rMVel = x2x2+x2y2+x2z2+y2x2+y2y2+y2z2+z2x2+z2y2+z2z2

          rFinal(iZeta,iComp,ipa,ipb) = Fact*rMVel
        end do

      end do
    end do
  end do
end do

return

end subroutine CmbnMV
