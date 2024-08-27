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
subroutine CmbnAC(Rnxyz,nZeta,la,lb,rKappa,final,Alpha,IfGrad,ld,nVecAC)
!***********************************************************************
!                                                                      *
! Object: compute the gradient of the overlap matrix.                  *
!                                                                      *
!     Author: Roland Lindh, Dept. of Theoretical Chemistry,            *
!             University of Lund, SWEDEN                               *
!             October '91.                                             *
!***********************************************************************

use Constants, only: Two

implicit none
integer nZeta, la, lb, ld, nVecAC
real*8 final(nZeta,(la+1)*(la+2)/2,(lb+1)*(lb+2)/2,4), rKappa(nZeta), Rnxyz(nZeta,3,0:la+ld,0:lb), Alpha(nZeta)
logical IfGrad(3)
integer ixa, ixb, iya, iyb, iza, izb, iZeta, iyaMax, iybMax, ipa, ipb
integer ixyz, ix, iz, Ind
real*8 tTwo, XA, YA, ZA
! Statement function for Cartesian index
Ind(ixyz,ix,iz) = (ixyz-ix)*(ixyz-ix+1)/2+iz+1

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
      ipa = Ind(la,ixa,iza)
      do iyb=0,iybMax
        izb = lb-ixb-iyb
        ipb = Ind(lb,ixb,izb)

        ! Combine overlap integral gradients

        nVecAC = 1
        do iZeta=1,nZeta
          final(iZeta,ipa,ipb,nVecAC) = rKappa(iZeta)*Rnxyz(iZeta,1,ixa,ixb)*Rnxyz(iZeta,2,iya,iyb)*Rnxyz(iZeta,3,iza,izb)
        end do
        tTwo = Two
        if (IfGrad(1)) then
          nVecAC = nVecAC+1
          if (ixa > 0) then
            xa = dble(-ixa)
            do iZeta=1,nZeta
              final(iZeta,ipa,ipb,nVecAC) = rKappa(iZeta)* &
                                            (tTwo*Alpha(iZeta)*Rnxyz(iZeta,1,ixa+1,ixb)+xa*Rnxyz(iZeta,1,ixa-1,ixb))* &
                                            Rnxyz(iZeta,2,iya,iyb)*Rnxyz(iZeta,3,iza,izb)
            end do
          else
            do iZeta=1,nZeta
              final(iZeta,ipa,ipb,nVecAC) = rKappa(iZeta)*tTwo*Alpha(iZeta)*Rnxyz(iZeta,1,ixa+1,ixb)*Rnxyz(iZeta,2,iya,iyb)* &
                                            Rnxyz(iZeta,3,iza,izb)
            end do
          end if
        end if
        if (IfGrad(2)) then
          nVecAC = nVecAC+1
          if (iya > 0) then
            ya = dble(-iya)
            do iZeta=1,nZeta
              final(iZeta,ipa,ipb,nVecAC) = rKappa(iZeta)*Rnxyz(iZeta,1,ixa,ixb)* &
                                            (tTwo*Alpha(iZeta)*Rnxyz(iZeta,2,iya+1,iyb)+ya*Rnxyz(iZeta,2,iya-1,iyb))* &
                                            Rnxyz(iZeta,3,iza,izb)
            end do
          else
            do iZeta=1,nZeta
              final(iZeta,ipa,ipb,nVecAC) = rKappa(iZeta)*Rnxyz(iZeta,1,ixa,ixb)*tTwo*Alpha(iZeta)*Rnxyz(iZeta,2,iya+1,iyb)* &
                                            Rnxyz(iZeta,3,iza,izb)
            end do
          end if
        end if
        if (IfGrad(3)) then
          nVecAC = nVecAC+1
          if (iza > 0) then
            za = dble(-iza)
            do iZeta=1,nZeta
              final(iZeta,ipa,ipb,nVecAC) = rKappa(iZeta)*Rnxyz(iZeta,1,ixa,ixb)*Rnxyz(iZeta,2,iya,iyb)* &
                                            (tTwo*Alpha(iZeta)*Rnxyz(iZeta,3,iza+1,izb)+za*Rnxyz(iZeta,3,iza-1,izb))
            end do
          else
            do iZeta=1,nZeta
              final(iZeta,ipa,ipb,nVecAC) = rKappa(iZeta)*Rnxyz(iZeta,1,ixa,ixb)*Rnxyz(iZeta,2,iya,iyb)*Two*Alpha(iZeta)* &
                                            Rnxyz(iZeta,3,iza+1,izb)
            end do
          end if
        end if

      end do
    end do
  end do
end do

return

end subroutine CmbnAC
