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

subroutine CmbnKE(Rnxyz,nZeta,la,lb,lr,Zeta,rKappa,rFinal,nComp,Txyz)
!***********************************************************************
!     Author: Roland Lindh, Dept. of Theoretical Chemistry,            *
!             University of Lund, SWEDEN                               *
!***********************************************************************

use Index_Functions, only: C_Ind, nTri_Elem1
use Constants, only: OneHalf
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: nZeta, la, lb, lr, nComp
real(kind=wp), intent(in) :: Rnxyz(nZeta,3,0:la+1,0:lb+1,0:lr), Zeta(nZeta), rKappa(nZeta), Txyz(nZeta,3,0:la,0:lb)
real(kind=wp), intent(out) :: rFinal(nZeta,nComp,nTri_Elem1(la),nTri_Elem1(lb))
integer(kind=iwp) :: iComp, ipa, ipb, ixa, ixb, iya, iyaMax, iyb, iybMax, iza, izb

!iRout = 134
!iPrint = nPrint(iRout)

iComp = 1
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

        ! Combine integrals

        rFinal(:,iComp,ipa,ipb) = rKappa*Zeta**(-OneHalf)*(Txyz(:,1,ixa,ixb)*Rnxyz(:,2,iya,iyb,0)*Rnxyz(:,3,iza,izb,0)+ &
                                                           Rnxyz(:,1,ixa,ixb,0)*Txyz(:,2,iya,iyb)*Rnxyz(:,3,iza,izb,0)+ &
                                                           Rnxyz(:,1,ixa,ixb,0)*Rnxyz(:,2,iya,iyb,0)*Txyz(:,3,iza,izb))

      end do
    end do
  end do
end do

return

end subroutine CmbnKE
