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
!               1991, Anders Bernhardsson                              *
!***********************************************************************

subroutine CmbnT1_mck(Rnxyz,nZeta,la,lb,Zeta,rKappa,rFinal,Txyz,Alpha,Beta,IfGrad)
!***********************************************************************
!     Author: Roland Lindh, Dept. of Theoretical Chemistry,            *
!             University of Lund, SWEDEN                               *
!             October '91                                              *
!             Anders Bernhardsson, Dept. of Theoretical Chemistry,     *
!             University of Lund, SWEDEN                               *
!             October '91                                              *
!***********************************************************************

use Index_Functions, only: C_Ind, nTri_Elem1
use Constants, only: Two, OneHalf
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp) :: nZeta, la, lb
real(kind=wp) :: Rnxyz(nZeta,3,0:la+2,0:lb+2), Zeta(nZeta), rKappa(nZeta), rFinal(nZeta,nTri_Elem1(la),nTri_Elem1(lb),1), &
                 Txyz(nZeta,3,0:la+1,0:lb+1), Alpha(nZeta), Beta(nZeta)
logical(kind=iwp) :: IfGrad(3,2)
integer(kind=iwp) :: ipa, ipb, ixa, ixb, iya, iyaMax, iyb, iybMax, iza, izb, iZeta
real(kind=wp) :: xa, xb, ya, yb, za, zb

!iRout = 134
!iPrint = nPrint(iRout)

!ii = nTri3_Elem(la)
!jj = nTri3_Elem(lb)
do iZeta=1,nZeta
  rKappa(iZeta) = rKappa(iZeta)*Zeta(iZeta)**(-OneHalf)
end do
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

        ! Combine integrals

        if (IfGrad(1,1)) then
          if (ixa > 0) then
            xa = real(-ixa,kind=wp)
            do iZeta=1,nZeta
              rFinal(iZeta,ipa,ipb,1) = rKappa(iZeta)*((Two*Txyz(iZeta,1,ixa+1,ixb)*Alpha(iZeta)+xa*Txyz(iZeta,1,ixa-1,ixb))* &
                                                       Rnxyz(iZeta,2,iya,iyb)*Rnxyz(iZeta,3,iza,izb)+ &
                                                       (Two*Rnxyz(iZeta,1,ixa+1,ixb)*Alpha(iZeta)+xa*Rnxyz(iZeta,1,ixa-1,ixb))* &
                                                       Txyz(iZeta,2,iya,iyb)*Rnxyz(iZeta,3,iza,izb)+ &
                                                       (Two*Rnxyz(iZeta,1,ixa+1,ixb)*Alpha(iZeta)+xa*Rnxyz(iZeta,1,ixa-1,ixb))* &
                                                       Rnxyz(iZeta,2,iya,iyb)*Txyz(iZeta,3,iza,izb))
            end do
          else
            do iZeta=1,nZeta
              rFinal(iZeta,ipa,ipb,1) = rKappa(iZeta)*Alpha(iZeta)*(Two*Txyz(iZeta,1,ixa+1,ixb)*Rnxyz(iZeta,2,iya,iyb)* &
                                                                    Rnxyz(iZeta,3,iza,izb)+Two*Rnxyz(iZeta,1,ixa+1,ixb)* &
                                                                    Txyz(iZeta,2,iya,iyb)*Rnxyz(iZeta,3,iza,izb)+ &
                                                                    Two*Rnxyz(iZeta,1,ixa+1,ixb)*Rnxyz(iZeta,2,iya,iyb)* &
                                                                    Txyz(iZeta,3,iza,izb))
            end do
          end if
        end if
        if (IfGrad(1,2)) then
          if (ixb > 0) then
            xb = real(-ixb,kind=wp)
            do iZeta=1,nZeta
              rFinal(iZeta,ipa,ipb,1) = rKappa(iZeta)*((Two*Txyz(iZeta,1,ixa,ixb+1)*Beta(iZeta)+xb*Txyz(iZeta,1,ixa,ixb-1))* &
                                                       Rnxyz(iZeta,2,iya,iyb)*Rnxyz(iZeta,3,iza,izb)+ &
                                                       (Two*Rnxyz(iZeta,1,ixa,ixb+1)*Beta(iZeta)+xb*Rnxyz(iZeta,1,ixa,ixb-1))* &
                                                       Txyz(iZeta,2,iya,iyb)*Rnxyz(iZeta,3,iza,izb)+ &
                                                       (Two*Rnxyz(iZeta,1,ixa,ixb+1)*Beta(iZeta)+xb*Rnxyz(iZeta,1,ixa,ixb-1))* &
                                                       Rnxyz(iZeta,2,iya,iyb)*Txyz(iZeta,3,iza,izb))
            end do
          else
            do iZeta=1,nZeta
              rFinal(iZeta,ipa,ipb,1) = rKappa(iZeta)*Beta(iZeta)*(Two*Txyz(iZeta,1,ixa,ixb+1)*Rnxyz(iZeta,2,iya,iyb)* &
                                                                   Rnxyz(iZeta,3,iza,izb)+Two*Rnxyz(iZeta,1,ixa,ixb+1)* &
                                                                   Txyz(iZeta,2,iya,iyb)*Rnxyz(iZeta,3,iza,izb)+ &
                                                                   Two*Rnxyz(iZeta,1,ixa,ixb+1)*Rnxyz(iZeta,2,iya,iyb)* &
                                                                   Txyz(iZeta,3,iza,izb))
            end do
          end if
        end if
        if (IfGrad(2,1)) then
          if (iya > 0) then
            ya = real(-iya,kind=wp)
            do iZeta=1,nZeta
              rFinal(iZeta,ipa,ipb,1) = rKappa(iZeta)*(Txyz(iZeta,1,ixa,ixb)*(Two*Rnxyz(iZeta,2,iya+1,iyb)*Alpha(iZeta)+ &
                                                       ya*Rnxyz(iZeta,2,iya-1,iyb))*Rnxyz(iZeta,3,iza,izb)+Rnxyz(iZeta,1,ixa,ixb)* &
                                                       (Two*Txyz(iZeta,2,iya+1,iyb)*Alpha(iZeta)+ya*Txyz(iZeta,2,iya-1,iyb))* &
                                                       Rnxyz(iZeta,3,iza,izb)+Rnxyz(iZeta,1,ixa,ixb)* &
                                                       (Two*Rnxyz(iZeta,2,iya+1,iyb)*Alpha(iZeta)+ya*Rnxyz(iZeta,2,iya-1,iyb))* &
                                                       Txyz(iZeta,3,iza,izb))
            end do
          else
            do iZeta=1,nZeta
              rFinal(iZeta,ipa,ipb,1) = rKappa(iZeta)*Alpha(iZeta)*(Txyz(iZeta,1,ixa,ixb)*Two*Rnxyz(iZeta,2,iya+1,iyb)* &
                                                                    Rnxyz(iZeta,3,iza,izb)+Rnxyz(iZeta,1,ixa,ixb)* &
                                                                    Two*Txyz(iZeta,2,iya+1,iyb)*Rnxyz(iZeta,3,iza,izb)+ &
                                                                    Rnxyz(iZeta,1,ixa,ixb)*Two*Rnxyz(iZeta,2,iya+1,iyb)* &
                                                                    Txyz(iZeta,3,iza,izb))
            end do
          end if
        end if
        if (IfGrad(2,2)) then
          if (iyb > 0) then
            yb = real(-iyb,kind=wp)
            do iZeta=1,nZeta
              rFinal(iZeta,ipa,ipb,1) = rKappa(iZeta)*(Txyz(iZeta,1,ixa,ixb)*(Two*Rnxyz(iZeta,2,iya,iyb+1)*Beta(iZeta)+ &
                                               yb*Rnxyz(iZeta,2,iya,iyb-1))*Rnxyz(iZeta,3,iza,izb)+Rnxyz(iZeta,1,ixa,ixb)* &
                                               (Two*Txyz(iZeta,2,iya,iyb+1)*Beta(iZeta)+yb*Txyz(iZeta,2,iya,iyb-1))* &
                                               Rnxyz(iZeta,3,iza,izb)+Rnxyz(iZeta,1,ixa,ixb)* &
                                               (Two*Rnxyz(iZeta,2,iya,iyb+1)*Beta(iZeta)+yb*Rnxyz(iZeta,2,iya,iyb-1))* &
                                               Txyz(iZeta,3,iza,izb))
            end do
          else
            do iZeta=1,nZeta
              rFinal(iZeta,ipa,ipb,1) = rKappa(iZeta)*Beta(iZeta)*(Txyz(iZeta,1,ixa,ixb)*Two*Rnxyz(iZeta,2,iya,iyb+1)* &
                                                                   Rnxyz(iZeta,3,iza,izb)+Rnxyz(iZeta,1,ixa,ixb)* &
                                                                   Two*Txyz(iZeta,2,iya,iyb+1)*Rnxyz(iZeta,3,iza,izb)+ &
                                                                   Rnxyz(iZeta,1,ixa,ixb)*Two*Rnxyz(iZeta,2,iya,iyb+1)* &
                                                                   Txyz(iZeta,3,iza,izb))
            end do
          end if
        end if
        if (IfGrad(3,1)) then
          if (iza > 0) then
            za = real(-iza,kind=wp)
            do iZeta=1,nZeta
              rFinal(iZeta,ipa,ipb,1) = rKappa(iZeta)*(Txyz(iZeta,1,ixa,ixb)*Rnxyz(iZeta,2,iya,iyb)* &
                                                       (Two*Rnxyz(iZeta,3,iza+1,izb)*Alpha(iZeta)+za*Rnxyz(iZeta,3,iza-1,izb))+ &
                                                       Rnxyz(iZeta,1,ixa,ixb)*Txyz(iZeta,2,iya,iyb)* &
                                                       (Two*Rnxyz(iZeta,3,iza+1,izb)*Alpha(iZeta)+za*Rnxyz(iZeta,3,iza-1,izb))+ &
                                                       Rnxyz(iZeta,1,ixa,ixb)*Rnxyz(iZeta,2,iya,iyb)* &
                                                       (Two*Txyz(iZeta,3,iza+1,izb)*Alpha(iZeta)+za*Txyz(iZeta,3,iza-1,izb)))
            end do
          else
            do iZeta=1,nZeta
              rFinal(iZeta,ipa,ipb,1) = rKappa(iZeta)*Alpha(iZeta)*(Txyz(iZeta,1,ixa,ixb)*Rnxyz(iZeta,2,iya,iyb)* &
                                                                    Two*Rnxyz(iZeta,3,iza+1,izb)+Rnxyz(iZeta,1,ixa,ixb)* &
                                                                    Txyz(iZeta,2,iya,iyb)*Two*Rnxyz(iZeta,3,iza+1,izb)+ &
                                                                    Rnxyz(iZeta,1,ixa,ixb)*Rnxyz(iZeta,2,iya,iyb)* &
                                                                    Two*Txyz(iZeta,3,iza+1,izb))
            end do
          end if
        end if
        if (IfGrad(3,2)) then
          if (izb > 0) then
            zb = real(-izb,kind=wp)
            do iZeta=1,nZeta
              rFinal(iZeta,ipa,ipb,1) = rKappa(iZeta)*(Txyz(iZeta,1,ixa,ixb)*Rnxyz(iZeta,2,iya,iyb)* &
                                                       (Two*Rnxyz(iZeta,3,iza,izb+1)*Beta(iZeta)+zb*Rnxyz(iZeta,3,iza,izb-1))+ &
                                                       Rnxyz(iZeta,1,ixa,ixb)*Txyz(iZeta,2,iya,iyb)* &
                                                       (Two*Rnxyz(iZeta,3,iza,izb+1)*Beta(iZeta)+zb*Rnxyz(iZeta,3,iza,izb-1))+ &
                                                       Rnxyz(iZeta,1,ixa,ixb)*Rnxyz(iZeta,2,iya,iyb)* &
                                                       (Two*Txyz(iZeta,3,iza,izb+1)*Beta(iZeta)+zb*Txyz(iZeta,3,iza,izb-1)))
            end do
          else
            do iZeta=1,nZeta
              rFinal(iZeta,ipa,ipb,1) = rKappa(iZeta)*Beta(iZeta)*(Txyz(iZeta,1,ixa,ixb)*Rnxyz(iZeta,2,iya,iyb)* &
                                                                   Two*Rnxyz(iZeta,3,iza,izb+1)+Rnxyz(iZeta,1,ixa,ixb)* &
                                                                   Txyz(iZeta,2,iya,iyb)*Two*Rnxyz(iZeta,3,iza,izb+1)+ &
                                                                   Rnxyz(iZeta,1,ixa,ixb)*Rnxyz(iZeta,2,iya,iyb)* &
                                                                   Two*Txyz(iZeta,3,iza,izb+1))
            end do
          end if
        end if

      end do
    end do
  end do
end do

return

end subroutine CmbnT1_mck
