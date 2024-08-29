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

use Constants, only: Two
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: nZeta, la, lb, ld
real(kind=wp), intent(in) :: Rnxyz(nZeta,3,0:la,0:lb+ld), rKappa(nZeta), Beta(nZeta)
real(kind=wp), intent(inout) :: rFinal(nZeta,(la+1)*(la+2)/2,(lb+1)*(lb+2)/2,4)
logical(kind=iwp), intent(in) :: IfGrad(3)
integer(kind=iwp), intent(out) :: nVecCB
integer(kind=iwp) :: ipa, ipb, ixa, ixb, iya, iyaMax, iyb, iybMax, iza, izb, iZeta
real(kind=wp) :: xB, yB, zB
! Statement function for Cartesian index
integer(kind=iwp) :: ixyz, ix, iz, Ind
Ind(ixyz,ix,iz) = (ixyz-ix)*(ixyz-ix+1)/2+iz+1

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
      ipa = Ind(la,ixa,iza)
      do iyb=0,iybMax
        izb = lb-ixb-iyb
        ipb = Ind(lb,ixb,izb)

        ! Combine overlap integral gradients

        nVecCB = 1
        do iZeta=1,nZeta
          rFinal(iZeta,ipa,ipb,nVecCB) = rKappa(iZeta)*Rnxyz(iZeta,1,ixa,ixb)*Rnxyz(iZeta,2,iya,iyb)*Rnxyz(iZeta,3,iza,izb)
        end do
        if (IfGrad(1)) then
          nVecCB = nVecCB+1
          if (ixb > 0) then
            xb = real(-ixb,kind=wp)
            do iZeta=1,nZeta
              rFinal(iZeta,ipa,ipb,nVecCB) = rKappa(iZeta)*(Two*Beta(iZeta)*Rnxyz(iZeta,1,ixa,ixb+1)+xb*Rnxyz(iZeta,1,ixa,ixb-1))* &
                                             Rnxyz(iZeta,2,iya,iyb)*Rnxyz(iZeta,3,iza,izb)
            end do
          else
            do iZeta=1,nZeta
              rFinal(iZeta,ipa,ipb,nVecCB) = rKappa(iZeta)*Two*Beta(iZeta)*Rnxyz(iZeta,1,ixa,ixb+1)*Rnxyz(iZeta,2,iya,iyb)* &
                                             Rnxyz(iZeta,3,iza,izb)
            end do
          end if
        end if
        if (IfGrad(2)) then
          nVecCB = nVecCB+1
          if (iyb > 0) then
            yb = real(-iyb,kind=wp)
            do iZeta=1,nZeta
              rFinal(iZeta,ipa,ipb,nVecCB) = rKappa(iZeta)*Rnxyz(iZeta,1,ixa,ixb)* &
                                             (Two*Beta(iZeta)*Rnxyz(iZeta,2,iya,iyb+1)+yb*Rnxyz(iZeta,2,iya,iyb-1))* &
                                             Rnxyz(iZeta,3,iza,izb)
            end do
          else
            do iZeta=1,nZeta
              rFinal(iZeta,ipa,ipb,nVecCB) = rKappa(iZeta)*Rnxyz(iZeta,1,ixa,ixb)*Two*Beta(iZeta)*Rnxyz(iZeta,2,iya,iyb+1)* &
                                             Rnxyz(iZeta,3,iza,izb)
            end do
          end if
        end if
        if (IfGrad(3)) then
          nVecCB = nVecCB+1
          if (izb > 0) then
            zb = real(-izb,kind=wp)
            do iZeta=1,nZeta
              rFinal(iZeta,ipa,ipb,nVecCB) = rKappa(iZeta)*Rnxyz(iZeta,1,ixa,ixb)*Rnxyz(iZeta,2,iya,iyb)*(Two*Beta(iZeta)* &
                                             Rnxyz(iZeta,3,iza,izb+1)+zb*Rnxyz(iZeta,3,iza,izb-1))
            end do
          else
            do iZeta=1,nZeta
              rFinal(iZeta,ipa,ipb,nVecCB) = rKappa(iZeta)*Rnxyz(iZeta,1,ixa,ixb)*Rnxyz(iZeta,2,iya,iyb)*Two*Beta(iZeta)* &
                                             Rnxyz(iZeta,3,iza,izb+1)
            end do
          end if
        end if

      end do
    end do
  end do
end do

return

end subroutine CmbnCB
