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

subroutine CmbnM2(Rnxyz,nZeta,la,lb,Zeta,rKappa,rFinal,Alpha,Beta,IfGrad,Fact,mVec)
!***********************************************************************
!                                                                      *
! Object: compute the gradient of the overlap matrix.                  *
!                                                                      *
!     Author: Roland Lindh, Dept. of Theoretical Chemistry,            *
!             University of Lund, SWEDEN                               *
!             October '91.                                             *
!***********************************************************************

use Index_Functions, only: C_Ind
use Constants, only: Two, Three
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: nZeta, la, lb
integer(kind=iwp), intent(out) :: mVec
logical(kind=iwp), intent(in) :: IfGrad(3,2)
real(kind=wp), intent(in) :: Rnxyz(nZeta,3,0:la+1,0:lb+1), Zeta(nZeta), Alpha(nZeta), Beta(nZeta), Fact
real(kind=wp), intent(inout) :: rKappa(nZeta)
real(kind=wp), intent(out) :: rFinal(nZeta,(la+1)*(la+2)/2,(lb+1)*(lb+2)/2,6)
integer(kind=iwp) :: ipa, ipb, iPrint, iRout, ixa, ixb, iya, iyaMax, iyb, iybMax, iza, izb, iZeta
real(kind=wp) :: xa, xb, ya, yb, za, zb
real(kind=wp), parameter :: exp32 = -Three/Two
#include "print.fh"

iRout = 134
iPrint = nPrint(iRout)

!ii = la*(la+1)*(la+2)/6
!jj = lb*(lb+1)*(lb+2)/6
do iZeta=1,nZeta
  rKappa(iZeta) = Fact*rKappa(iZeta)*Zeta(iZeta)**exp32
end do
if (iPrint >= 99) then
  call RecPrt(' In CmbnM2: Zeta  ',' ',Zeta,1,nZeta)
  call RecPrt(' In CmbnM2: rKappa',' ',rKappa,1,nZeta)
  call RecPrt(' In CmbnM2: Alpha ',' ',Alpha,1,nZeta)
  call RecPrt(' In CmbnM2: Beta  ',' ',Beta,1,nZeta)
  call RecPrt(' In CmbnM2: Rnxyz ',' ',Rnxyz,nZeta*3,(la+2)*(lb+2))
end if
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

        ! Combine overlap integrals

        mVec = 0
        if (IfGrad(1,1)) then
          mVec = mVec+1
          if (ixa > 0) then
            xa = -real(ixa,kind=wp)
            do iZeta=1,nZeta
              rFinal(iZeta,ipa,ipb,mVec) = rKappa(iZeta)* &
                         (Two*Alpha(iZeta)*Rnxyz(iZeta,1,ixa+1,ixb)+ &
                                        xa*Rnxyz(iZeta,1,ixa-1,ixb))* &
                                           Rnxyz(iZeta,2,iya,iyb)* &
                                           Rnxyz(iZeta,3,iza,izb)
            end do
          else
            do iZeta=1,nZeta
              rFinal(iZeta,ipa,ipb,mVec) = rKappa(iZeta)* &
                          Two*Alpha(iZeta)*Rnxyz(iZeta,1,ixa+1,ixb)* &
                                           Rnxyz(iZeta,2,iya,iyb)* &
                                           Rnxyz(iZeta,3,iza,izb)
            end do
          end if
        end if
        if (IfGrad(1,2)) then
          mVec = mVec+1
          if (ixb > 0) then
            xb = -real(ixb,kind=wp)
            do iZeta=1,nZeta
              rFinal(iZeta,ipa,ipb,mVec) = rKappa(iZeta)* &
                          (Two*Beta(iZeta)*Rnxyz(iZeta,1,ixa,ixb+1)+ &
                                        xb*Rnxyz(iZeta,1,ixa,ixb-1))* &
                                           Rnxyz(iZeta,2,iya,iyb)* &
                                           Rnxyz(iZeta,3,iza,izb)
            end do
          else
            do iZeta=1,nZeta
              rFinal(iZeta,ipa,ipb,mVec) = rKappa(iZeta)* &
                           Two*Beta(iZeta)*Rnxyz(iZeta,1,ixa,ixb+1)* &
                                           Rnxyz(iZeta,2,iya,iyb)* &
                                           Rnxyz(iZeta,3,iza,izb)
            end do
          end if
        end if
        if (IfGrad(2,1)) then
          mVec = mVec+1
          if (iya > 0) then
            ya = -real(iya,kind=wp)
            do iZeta=1,nZeta
              rFinal(iZeta,ipa,ipb,mVec) = rKappa(iZeta)* &
                                           Rnxyz(iZeta,1,ixa,ixb)* &
                         (Two*Alpha(iZeta)*Rnxyz(iZeta,2,iya+1,iyb)+ &
                                        ya*Rnxyz(iZeta,2,iya-1,iyb))* &
                                           Rnxyz(iZeta,3,iza,izb)
            end do
          else
            do iZeta=1,nZeta
              rFinal(iZeta,ipa,ipb,mVec) = rKappa(iZeta)* &
                                           Rnxyz(iZeta,1,ixa,ixb)* &
                          Two*Alpha(iZeta)*Rnxyz(iZeta,2,iya+1,iyb)* &
                                           Rnxyz(iZeta,3,iza,izb)
            end do
          end if
        end if
        if (IfGrad(2,2)) then
          mVec = mVec+1
          if (iyb > 0) then
            yb = -real(iyb,kind=wp)
            do iZeta=1,nZeta
              rFinal(iZeta,ipa,ipb,mVec) = rKappa(iZeta)* &
                                           Rnxyz(iZeta,1,ixa,ixb)* &
                          (Two*Beta(iZeta)*Rnxyz(iZeta,2,iya,iyb+1)+ &
                                        yb*Rnxyz(iZeta,2,iya,iyb-1))* &
                                           Rnxyz(iZeta,3,iza,izb)
            end do
          else
            do iZeta=1,nZeta
              rFinal(iZeta,ipa,ipb,mVec) = rKappa(iZeta)* &
                                           Rnxyz(iZeta,1,ixa,ixb)* &
                           Two*Beta(iZeta)*Rnxyz(iZeta,2,iya,iyb+1)* &
                                           Rnxyz(iZeta,3,iza,izb)
            end do
          end if
        end if
        if (IfGrad(3,1)) then
          mVec = mVec+1
          if (iza > 0) then
            za = -real(iza,kind=wp)
            do iZeta=1,nZeta
              rFinal(iZeta,ipa,ipb,mVec) = rKappa(iZeta)* &
                                           Rnxyz(iZeta,1,ixa,ixb)* &
                                           Rnxyz(iZeta,2,iya,iyb)* &
                         (Two*Alpha(iZeta)*Rnxyz(iZeta,3,iza+1,izb)+ &
                                        za*Rnxyz(iZeta,3,iza-1,izb))
            end do
          else
            do iZeta=1,nZeta
              rFinal(iZeta,ipa,ipb,mVec) = rKappa(iZeta)* &
                                           Rnxyz(iZeta,1,ixa,ixb)* &
                                           Rnxyz(iZeta,2,iya,iyb)* &
                          Two*Alpha(iZeta)*Rnxyz(iZeta,3,iza+1,izb)
            end do
          end if
        end if
        if (IfGrad(3,2)) then
          mVec = mVec+1
          if (izb > 0) then
            zb = -real(izb,kind=wp)
            do iZeta=1,nZeta
              rFinal(iZeta,ipa,ipb,mVec) = rKappa(iZeta)* &
                                           Rnxyz(iZeta,1,ixa,ixb)* &
                                           Rnxyz(iZeta,2,iya,iyb)* &
                          (Two*Beta(iZeta)*Rnxyz(iZeta,3,iza,izb+1)+ &
                                        zb*Rnxyz(iZeta,3,iza,izb-1))
            end do
          else
            do iZeta=1,nZeta
              rFinal(iZeta,ipa,ipb,mVec) = rKappa(iZeta)* &
                                           Rnxyz(iZeta,1,ixa,ixb)* &
                                           Rnxyz(iZeta,2,iya,iyb)* &
                           Two*Beta(iZeta)*Rnxyz(iZeta,3,iza,izb+1)
            end do
          end if
        end if

      end do
    end do
  end do
end do

return

end subroutine CmbnM2
