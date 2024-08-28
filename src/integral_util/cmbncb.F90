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
subroutine CmbnCB(Rnxyz,nZeta,la,lb,rKappa,final,Beta,IfGrad,ld,nVecCB)
!***********************************************************************
!                                                                      *
! Object: compute the gradient of the overlap matrix.                  *
!                                                                      *
!     Author: Roland Lindh, Dept. of Theoretical Chemistry,            *
!             University of Lund, SWEDEN                               *
!             October '91.                                             *
!***********************************************************************

use Constants, only: Two
use Definitions, only: wp

implicit none
integer nZeta, la, lb, ld, nVecCB
real*8 final(nZeta,(la+1)*(la+2)/2,(lb+1)*(lb+2)/2,4), rKappa(nZeta), Beta(nZeta), Rnxyz(nZeta,3,0:la,0:lb+ld)
logical IfGrad(3)
integer ixyz, ix, iz, Ind
integer ixa, ixb, iya, iyb, iza, izb, ipa, ipb, iZeta, iyaMax, iybMax
real*8 xB, yB, zB, tTwo
! Statement function for Cartesian index
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
          final(iZeta,ipa,ipb,nVecCB) = rKappa(iZeta)*Rnxyz(iZeta,1,ixa,ixb)*Rnxyz(iZeta,2,iya,iyb)*Rnxyz(iZeta,3,iza,izb)
        end do
        tTwo = Two
        if (IfGrad(1)) then
          nVecCB = nVecCB+1
          if (ixb > 0) then
            xb = real(-ixb,kind=wp)
            do iZeta=1,nZeta
              final(iZeta,ipa,ipb,nVecCB) = rKappa(iZeta)*(tTwo*Beta(iZeta)*Rnxyz(iZeta,1,ixa,ixb+1)+xb*Rnxyz(iZeta,1,ixa,ixb-1))* &
                                            Rnxyz(iZeta,2,iya,iyb)*Rnxyz(iZeta,3,iza,izb)
            end do
          else
            do iZeta=1,nZeta
              final(iZeta,ipa,ipb,nVecCB) = rKappa(iZeta)*tTwo*Beta(iZeta)*Rnxyz(iZeta,1,ixa,ixb+1)*Rnxyz(iZeta,2,iya,iyb)* &
                                            Rnxyz(iZeta,3,iza,izb)
            end do
          end if
        end if
        if (IfGrad(2)) then
          nVecCB = nVecCB+1
          if (iyb > 0) then
            yb = real(-iyb,kind=wp)
            do iZeta=1,nZeta
              final(iZeta,ipa,ipb,nVecCB) = rKappa(iZeta)*Rnxyz(iZeta,1,ixa,ixb)* &
                                            (tTwo*Beta(iZeta)*Rnxyz(iZeta,2,iya,iyb+1)+yb*Rnxyz(iZeta,2,iya,iyb-1))* &
                                            Rnxyz(iZeta,3,iza,izb)
            end do
          else
            do iZeta=1,nZeta
              final(iZeta,ipa,ipb,nVecCB) = rKappa(iZeta)*Rnxyz(iZeta,1,ixa,ixb)*tTwo*Beta(iZeta)*Rnxyz(iZeta,2,iya,iyb+1)* &
                                            Rnxyz(iZeta,3,iza,izb)
            end do
          end if
        end if
        if (IfGrad(3)) then
          nVecCB = nVecCB+1
          if (izb > 0) then
            zb = real(-izb,kind=wp)
            do iZeta=1,nZeta
              final(iZeta,ipa,ipb,nVecCB) = rKappa(iZeta)*Rnxyz(iZeta,1,ixa,ixb)*Rnxyz(iZeta,2,iya,iyb)*(tTwo*Beta(iZeta)* &
                                            Rnxyz(iZeta,3,iza,izb+1)+zb*Rnxyz(iZeta,3,iza,izb-1))
            end do
          else
            do iZeta=1,nZeta
              final(iZeta,ipa,ipb,nVecCB) = rKappa(iZeta)*Rnxyz(iZeta,1,ixa,ixb)*Rnxyz(iZeta,2,iya,iyb)*tTwo*Beta(iZeta)* &
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
