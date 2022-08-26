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
!               1995, Anders Bernhardsson                              *
!***********************************************************************

subroutine CmbnS1_mck(Rnxyz,nZeta,la,lb,Zeta,rKappa,rFinal,Alpha,Beta,IfGrad)
!***********************************************************************
!                                                                      *
! Object: compute the gradient of the overlap matrix.                  *
!                                                                      *
!     Author: Roland Lindh,                                            *
!             Dept. of Theoretical Chemistry,                          *
!             University of Lund, SWEDEN                               *
!             October '91.                                             *
!             Anders Bernhardsson                                      *
!             Dept. of Theoretical Chemistry,                          *
!             University of Lund, SWEDEN                               *
!              95.                                                     *
!***********************************************************************

use Index_Functions, only: C_Ind
use Constants, only: Two, OneHalf
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp) :: nZeta, la, lb
real(kind=wp) :: Rnxyz(nZeta,3,0:la+1,0:lb+1), Zeta(nZeta), rKappa(nZeta), rFinal(nZeta,(la+1)*(la+2)/2,(lb+1)*(lb+2)/2,1), &
                 Alpha(nZeta), Beta(nZeta)
logical(kind=iwp) :: IfGrad(3,2)
integer(kind=iwp) :: ipa, ipb, ixa, ixb, iya, iyaMax, iyb, iybMax, iza, izb, iZeta
real(kind=wp) :: xa, xb, ya, yb, za, zb

!iRout = 134
!iPrint = nPrint(iRout)

!ii = la*(la+1)*(la+2)/6
!jj = lb*(lb+1)*(lb+2)/6
do iZeta=1,nZeta
  rKappa(iZeta) = rKappa(iZeta)*Zeta(iZeta)**(-OneHalf)
end do
!if (iPrint >= 99) then
!  call RecPrt(' In CmbnS1_mck: Zeta  ',' ',Zeta,1,nZeta)
!  call RecPrt(' In CmbnS1_mck: rKappa',' ',rKappa,1,nZeta)
!  call RecPrt(' In CmbnS1_mck: Alpha ',' ',Alpha,1,nZeta)
!  call RecPrt(' In CmbnS1_mck: Beta  ',' ',Beta,1,nZeta)
!end if
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

        !write(u6,*) ' papb=', papb
        if (IfGrad(1,1)) then
          if (ixa > 0) then
            xa = real(-ixa,kind=wp)
            do iZeta=1,nZeta
              !rFinal(iZeta,ipa,ipb,1) = papb*rKappa(iZeta)*
              rFinal(iZeta,ipa,ipb,1) = rKappa(iZeta)*(Two*Alpha(iZeta)*Rnxyz(iZeta,1,ixa+1,ixb)+xa*Rnxyz(iZeta,1,ixa-1,ixb))* &
                                        Rnxyz(iZeta,2,iya,iyb)*Rnxyz(iZeta,3,iza,izb)
            end do
          else
            do iZeta=1,nZeta
              !rFinal(iZeta,ipa,ipb,1) = papb*rKappa(iZeta)*
              rFinal(iZeta,ipa,ipb,1) = rKappa(iZeta)*Two*Alpha(iZeta)*Rnxyz(iZeta,1,ixa+1,ixb)*Rnxyz(iZeta,2,iya,iyb)* &
                                        Rnxyz(iZeta,3,iza,izb)
            end do
          end if
        end if
        if (IfGrad(1,2)) then
          if (ixb > 0) then
            xb = real(-ixb,kind=wp)
            do iZeta=1,nZeta
              !rFinal(iZeta,ipa,ipb,1) = papb*rKappa(iZeta)*
              rFinal(iZeta,ipa,ipb,1) = rKappa(iZeta)*(Two*Beta(iZeta)*Rnxyz(iZeta,1,ixa,ixb+1)+xb*Rnxyz(iZeta,1,ixa,ixb-1))* &
                                        Rnxyz(iZeta,2,iya,iyb)*Rnxyz(iZeta,3,iza,izb)
            end do
          else
            do iZeta=1,nZeta
              !rFinal(iZeta,ipa,ipb,1) = papb*rKappa(iZeta)*
              rFinal(iZeta,ipa,ipb,1) = rKappa(iZeta)*Two*Beta(iZeta)*Rnxyz(iZeta,1,ixa,ixb+1)*Rnxyz(iZeta,2,iya,iyb)* &
                                        Rnxyz(iZeta,3,iza,izb)
            end do
          end if
        end if
        if (IfGrad(2,1)) then
          if (iya > 0) then
            ya = real(-iya,kind=wp)
            do iZeta=1,nZeta
              !rFinal(iZeta,ipa,ipb,1) = papb*rKappa(iZeta)*
              rFinal(iZeta,ipa,ipb,1) = rKappa(iZeta)*Rnxyz(iZeta,1,ixa,ixb)*(Two*Alpha(iZeta)*Rnxyz(iZeta,2,iya+1,iyb)+ &
                                        ya*Rnxyz(iZeta,2,iya-1,iyb))*Rnxyz(iZeta,3,iza,izb)
            end do
          else
            do iZeta=1,nZeta
              !rFinal(iZeta,ipa,ipb,1) = papb*rKappa(iZeta)*
              rFinal(iZeta,ipa,ipb,1) = rKappa(iZeta)*Rnxyz(iZeta,1,ixa,ixb)*Two*Alpha(iZeta)*Rnxyz(iZeta,2,iya+1,iyb)* &
                                        Rnxyz(iZeta,3,iza,izb)
            end do
          end if
        end if
        if (IfGrad(2,2)) then
          if (iyb > 0) then
            yb = real(-iyb,kind=wp)
            do iZeta=1,nZeta
              !rFinal(iZeta,ipa,ipb,1) = papb*rKappa(iZeta)*
              rFinal(iZeta,ipa,ipb,1) = rKappa(iZeta)*Rnxyz(iZeta,1,ixa,ixb)*(Two*Beta(iZeta)*Rnxyz(iZeta,2,iya,iyb+1)+ &
                                        yb*Rnxyz(iZeta,2,iya,iyb-1))*Rnxyz(iZeta,3,iza,izb)
            end do
          else
            do iZeta=1,nZeta
              !rFinal(iZeta,ipa,ipb,1) = papb*rKappa(iZeta)*
              rFinal(iZeta,ipa,ipb,1) = rKappa(iZeta)*Rnxyz(iZeta,1,ixa,ixb)*Two*Beta(iZeta)*Rnxyz(iZeta,2,iya,iyb+1)* &
                                        Rnxyz(iZeta,3,iza,izb)
            end do
          end if
        end if
        if (IfGrad(3,1)) then
          if (iza > 0) then
            za = real(-iza,kind=wp)
            do iZeta=1,nZeta
              !rFinal(iZeta,ipa,ipb,1) = papb*rKappa(iZeta)*
              rFinal(iZeta,ipa,ipb,1) = rKappa(iZeta)*Rnxyz(iZeta,1,ixa,ixb)*Rnxyz(iZeta,2,iya,iyb)*(Two*Alpha(iZeta)* &
                                        Rnxyz(iZeta,3,iza+1,izb)+za*Rnxyz(iZeta,3,iza-1,izb))
            end do
          else
            do iZeta=1,nZeta
              !rFinal(iZeta,ipa,ipb,1) = papb*rKappa(iZeta)*
              rFinal(iZeta,ipa,ipb,1) = rKappa(iZeta)*Rnxyz(iZeta,1,ixa,ixb)*Rnxyz(iZeta,2,iya,iyb)*Two*Alpha(iZeta)* &
                                        Rnxyz(iZeta,3,iza+1,izb)
            end do
          end if
        end if
        if (IfGrad(3,2)) then
          if (izb > 0) then
            zb = real(-izb,kind=wp)
            do iZeta=1,nZeta
              !rFinal(iZeta,ipa,ipb,1) = papb*rKappa(iZeta)*
              rFinal(iZeta,ipa,ipb,1) = rKappa(iZeta)*Rnxyz(iZeta,1,ixa,ixb)*Rnxyz(iZeta,2,iya,iyb)*(Two*Beta(iZeta)* &
                                        Rnxyz(iZeta,3,iza,izb+1)+zb*Rnxyz(iZeta,3,iza,izb-1))
            end do
          else
            do iZeta=1,nZeta
              !rFinal(iZeta,ipa,ipb,1) = papb*rKappa(iZeta)*
              rFinal(iZeta,ipa,ipb,1) = rKappa(iZeta)*Rnxyz(iZeta,1,ixa,ixb)*Rnxyz(iZeta,2,iya,iyb)*Two*Beta(iZeta)* &
                                        Rnxyz(iZeta,3,iza,izb+1)
            end do
          end if
        end if

      end do
    end do
  end do
end do

return

end subroutine CmbnS1_mck
