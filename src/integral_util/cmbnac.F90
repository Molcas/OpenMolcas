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
subroutine CmbnAC(Rnxyz,nZeta,la,lb,rKappa,rFinal,Alpha,IfGrad,ld,nVecAC)
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
real(kind=wp), intent(in) :: Rnxyz(nZeta,3,0:la+ld,0:lb), rKappa(nZeta), Alpha(nZeta)
real(kind=wp), intent(inout) :: rFinal(nZeta,nTri_Elem1(la),nTri_Elem1(lb),4)
logical(kind=iwp), intent(in) :: IfGrad(3)
integer(kind=iwp), intent(out) :: nVecAC
integer(kind=iwp) :: ipa, ipb, ixa, ixb, iya, iyaMax, iyb, iybMax, iza, izb
real(kind=wp) :: xA, yA, zA

#ifdef _DEBUGPRINT_
call RecPrt(' In CmbnAC: rKappa',' ',rKappa,1,nZeta)
call RecPrt(' In CmbnAC: Alpha ',' ',Alpha,1,nZeta)
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

        nVecAC = 1
        rFinal(:,ipa,ipb,nVecAC) = rKappa(:)*Rnxyz(:,1,ixa,ixb)*Rnxyz(:,2,iya,iyb)*Rnxyz(:,3,iza,izb)
        if (IfGrad(1)) then
          nVecAC = nVecAC+1
          if (ixa > 0) then
            xa = real(-ixa,kind=wp)
            rFinal(:,ipa,ipb,nVecAC) = rKappa(:)*(Two*Alpha(:)*Rnxyz(:,1,ixa+1,ixb)+xa*Rnxyz(:,1,ixa-1,ixb))* &
                                       Rnxyz(:,2,iya,iyb)*Rnxyz(:,3,iza,izb)
          else
            rFinal(:,ipa,ipb,nVecAC) = rKappa(:)*Two*Alpha(:)*Rnxyz(:,1,ixa+1,ixb)*Rnxyz(:,2,iya,iyb)*Rnxyz(:,3,iza,izb)
          end if
        end if
        if (IfGrad(2)) then
          nVecAC = nVecAC+1
          if (iya > 0) then
            ya = real(-iya,kind=wp)
            rFinal(:,ipa,ipb,nVecAC) = rKappa(:)*Rnxyz(:,1,ixa,ixb)*(Two*Alpha(:)*Rnxyz(:,2,iya+1,iyb)+ya*Rnxyz(:,2,iya-1,iyb))* &
                                       Rnxyz(:,3,iza,izb)
          else
            rFinal(:,ipa,ipb,nVecAC) = rKappa(:)*Rnxyz(:,1,ixa,ixb)*Two*Alpha(:)*Rnxyz(:,2,iya+1,iyb)*Rnxyz(:,3,iza,izb)
          end if
        end if
        if (IfGrad(3)) then
          nVecAC = nVecAC+1
          if (iza > 0) then
            za = real(-iza,kind=wp)
            rFinal(:,ipa,ipb,nVecAC) = rKappa(:)*Rnxyz(:,1,ixa,ixb)*Rnxyz(:,2,iya,iyb)* &
                                       (Two*Alpha(:)*Rnxyz(:,3,iza+1,izb)+za*Rnxyz(:,3,iza-1,izb))
          else
            rFinal(:,ipa,ipb,nVecAC) = rKappa(:)*Rnxyz(:,1,ixa,ixb)*Rnxyz(:,2,iya,iyb)*Two*Alpha(:)*Rnxyz(:,3,iza+1,izb)
          end if
        end if

      end do
    end do
  end do
end do

return

end subroutine CmbnAC
