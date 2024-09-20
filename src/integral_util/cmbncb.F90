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
subroutine CmbnCB(Rnxyz,nZeta,la,lb,rKappa,rFinal,Beta,IfGrad,ld,nVecCB)
!***********************************************************************
!                                                                      *
! Object: compute the gradient of the overlap matrix.                  *
!                                                                      *
!     Author: Roland Lindh, Dept. of Theoretical Chemistry,            *
!             University of Lund, SWEDEN                               *
!             October '91.                                             *
!***********************************************************************

use Index_Functions, only: C_Ind, nTri_Elem1
use Constants, only: Two
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: nZeta, la, lb, ld
real(kind=wp), intent(in) :: Rnxyz(nZeta,3,0:la,0:lb+ld), rKappa(nZeta), Beta(nZeta)
real(kind=wp), intent(inout) :: rFinal(nZeta,nTri_Elem1(la),nTri_Elem1(lb),4)
logical(kind=iwp), intent(in) :: IfGrad(3)
integer(kind=iwp), intent(out) :: nVecCB
integer(kind=iwp) :: ipa, ipb, ixa, ixb, iya, iyaMax, iyb, iybMax, iza, izb
real(kind=wp) :: xB, yB, zB

#ifdef _DEBUGPRINT_
call RecPrt(' In CmbnCB: rKappa',' ',rKappa,1,nZeta)
call RecPrt(' In CmbnCB: Beta  ',' ',Beta,1,nZeta)
#endif
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

        ! Combine overlap integral gradients

        nVecCB = 1
        rFinal(:,ipa,ipb,nVecCB) = rKappa(:)*Rnxyz(:,1,ixa,ixb)*Rnxyz(:,2,iya,iyb)*Rnxyz(:,3,iza,izb)
        if (IfGrad(1)) then
          nVecCB = nVecCB+1
          if (ixb > 0) then
            xb = real(-ixb,kind=wp)
            rFinal(:,ipa,ipb,nVecCB) = rKappa(:)*(Two*Beta(:)*Rnxyz(:,1,ixa,ixb+1)+xb*Rnxyz(:,1,ixa,ixb-1))* &
                                       Rnxyz(:,2,iya,iyb)*Rnxyz(:,3,iza,izb)
          else
            rFinal(:,ipa,ipb,nVecCB) = rKappa(:)*Two*Beta(:)*Rnxyz(:,1,ixa,ixb+1)*Rnxyz(:,2,iya,iyb)*Rnxyz(:,3,iza,izb)
          end if
        end if
        if (IfGrad(2)) then
          nVecCB = nVecCB+1
          if (iyb > 0) then
            yb = real(-iyb,kind=wp)
            rFinal(:,ipa,ipb,nVecCB) = rKappa(:)*Rnxyz(:,1,ixa,ixb)*(Two*Beta(:)*Rnxyz(:,2,iya,iyb+1)+yb*Rnxyz(:,2,iya,iyb-1))* &
                                       Rnxyz(:,3,iza,izb)
          else
            rFinal(:,ipa,ipb,nVecCB) = rKappa(:)*Rnxyz(:,1,ixa,ixb)*Two*Beta(:)*Rnxyz(:,2,iya,iyb+1)*Rnxyz(:,3,iza,izb)
          end if
        end if
        if (IfGrad(3)) then
          nVecCB = nVecCB+1
          if (izb > 0) then
            zb = real(-izb,kind=wp)
            rFinal(:,ipa,ipb,nVecCB) = rKappa(:)*Rnxyz(:,1,ixa,ixb)*Rnxyz(:,2,iya,iyb)*(Two*Beta(:)* &
                                       Rnxyz(:,3,iza,izb+1)+zb*Rnxyz(:,3,iza,izb-1))
          else
            rFinal(:,ipa,ipb,nVecCB) = rKappa(:)*Rnxyz(:,1,ixa,ixb)*Rnxyz(:,2,iya,iyb)*Two*Beta(:)*Rnxyz(:,3,iza,izb+1)
          end if
        end if

      end do
    end do
  end do
end do

return

end subroutine CmbnCB
