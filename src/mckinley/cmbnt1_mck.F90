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
integer(kind=iwp), intent(in) :: nZeta, la, lb
real(kind=wp), intent(in) :: Rnxyz(nZeta,3,0:la+2,0:lb+2), Zeta(nZeta), Txyz(nZeta,3,0:la+1,0:lb+1), Alpha(nZeta), Beta(nZeta)
real(kind=wp), intent(inout) :: rKappa(nZeta), rFinal(nZeta,nTri_Elem1(la),nTri_Elem1(lb),1)
logical(kind=iwp), intent(in) :: IfGrad(3,2)
integer(kind=iwp) :: ipa, ipb, ixa, ixb, iya, iyaMax, iyb, iybMax, iza, izb
real(kind=wp) :: xa, xb, ya, yb, za, zb

!iRout = 134
!iPrint = nPrint(iRout)

!ii = nTri3_Elem(la)
!jj = nTri3_Elem(lb)
rKappa(:) = rKappa*Zeta**(-OneHalf)
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
            rFinal(:,ipa,ipb,1) = rKappa(:)*((Two*Txyz(:,1,ixa+1,ixb)*Alpha(:)+xa*Txyz(:,1,ixa-1,ixb))* &
                                             Rnxyz(:,2,iya,iyb)*Rnxyz(:,3,iza,izb)+ &
                                             (Two*Rnxyz(:,1,ixa+1,ixb)*Alpha(:)+xa*Rnxyz(:,1,ixa-1,ixb))* &
                                             Txyz(:,2,iya,iyb)*Rnxyz(:,3,iza,izb)+ &
                                             (Two*Rnxyz(:,1,ixa+1,ixb)*Alpha(:)+xa*Rnxyz(:,1,ixa-1,ixb))* &
                                             Rnxyz(:,2,iya,iyb)*Txyz(:,3,iza,izb))
          else
            rFinal(:,ipa,ipb,1) = rKappa(:)*Two*Alpha(:)*(Txyz(:,1,ixa+1,ixb)*Rnxyz(:,2,iya,iyb)*Rnxyz(:,3,iza,izb)+ &
                                                          Rnxyz(:,1,ixa+1,ixb)*Txyz(:,2,iya,iyb)*Rnxyz(:,3,iza,izb)+ &
                                                          Rnxyz(:,1,ixa+1,ixb)*Rnxyz(:,2,iya,iyb)*Txyz(:,3,iza,izb))
          end if
        end if
        if (IfGrad(1,2)) then
          if (ixb > 0) then
            xb = real(-ixb,kind=wp)
            rFinal(:,ipa,ipb,1) = rKappa(:)*((Two*Txyz(:,1,ixa,ixb+1)*Beta(:)+xb*Txyz(:,1,ixa,ixb-1))* &
                                             Rnxyz(:,2,iya,iyb)*Rnxyz(:,3,iza,izb)+ &
                                             (Two*Rnxyz(:,1,ixa,ixb+1)*Beta(:)+xb*Rnxyz(:,1,ixa,ixb-1))* &
                                             Txyz(:,2,iya,iyb)*Rnxyz(:,3,iza,izb)+ &
                                             (Two*Rnxyz(:,1,ixa,ixb+1)*Beta(:)+xb*Rnxyz(:,1,ixa,ixb-1))* &
                                             Rnxyz(:,2,iya,iyb)*Txyz(:,3,iza,izb))
          else
            rFinal(:,ipa,ipb,1) = rKappa(:)*Two*Beta(:)*(Txyz(:,1,ixa,ixb+1)*Rnxyz(:,2,iya,iyb)*Rnxyz(:,3,iza,izb)+ &
                                                         Rnxyz(:,1,ixa,ixb+1)*Txyz(:,2,iya,iyb)*Rnxyz(:,3,iza,izb)+ &
                                                         Rnxyz(:,1,ixa,ixb+1)*Rnxyz(:,2,iya,iyb)*Txyz(:,3,iza,izb))
          end if
        end if
        if (IfGrad(2,1)) then
          if (iya > 0) then
            ya = real(-iya,kind=wp)
            rFinal(:,ipa,ipb,1) = rKappa(:)*(Txyz(:,1,ixa,ixb)*(Two*Rnxyz(:,2,iya+1,iyb)*Alpha(:)+ya*Rnxyz(:,2,iya-1,iyb))* &
                                             Rnxyz(:,3,iza,izb)+ &
                                             Rnxyz(:,1,ixa,ixb)*(Two*Txyz(:,2,iya+1,iyb)*Alpha(:)+ya*Txyz(:,2,iya-1,iyb))* &
                                             Rnxyz(:,3,iza,izb)+ &
                                             Rnxyz(:,1,ixa,ixb)*(Two*Rnxyz(:,2,iya+1,iyb)*Alpha(:)+ya*Rnxyz(:,2,iya-1,iyb))* &
                                             Txyz(:,3,iza,izb))
          else
            rFinal(:,ipa,ipb,1) = rKappa(:)*Two*Alpha(:)*(Txyz(:,1,ixa,ixb)*Rnxyz(:,2,iya+1,iyb)*Rnxyz(:,3,iza,izb)+ &
                                                          Rnxyz(:,1,ixa,ixb)*Txyz(:,2,iya+1,iyb)*Rnxyz(:,3,iza,izb)+ &
                                                          Rnxyz(:,1,ixa,ixb)*Rnxyz(:,2,iya+1,iyb)*Txyz(:,3,iza,izb))
          end if
        end if
        if (IfGrad(2,2)) then
          if (iyb > 0) then
            yb = real(-iyb,kind=wp)
            rFinal(:,ipa,ipb,1) = rKappa(:)*(Txyz(:,1,ixa,ixb)*(Two*Rnxyz(:,2,iya,iyb+1)*Beta(:)+yb*Rnxyz(:,2,iya,iyb-1))* &
                                             Rnxyz(:,3,iza,izb)+ &
                                             Rnxyz(:,1,ixa,ixb)*(Two*Txyz(:,2,iya,iyb+1)*Beta(:)+yb*Txyz(:,2,iya,iyb-1))* &
                                             Rnxyz(:,3,iza,izb)+ &
                                             Rnxyz(:,1,ixa,ixb)*(Two*Rnxyz(:,2,iya,iyb+1)*Beta(:)+yb*Rnxyz(:,2,iya,iyb-1))* &
                                             Txyz(:,3,iza,izb))
          else
            rFinal(:,ipa,ipb,1) = rKappa(:)*Two*Beta(:)*(Txyz(:,1,ixa,ixb)*Rnxyz(:,2,iya,iyb+1)*Rnxyz(:,3,iza,izb)+ &
                                                         Rnxyz(:,1,ixa,ixb)*Txyz(:,2,iya,iyb+1)*Rnxyz(:,3,iza,izb)+ &
                                                         Rnxyz(:,1,ixa,ixb)*Rnxyz(:,2,iya,iyb+1)*Txyz(:,3,iza,izb))
          end if
        end if
        if (IfGrad(3,1)) then
          if (iza > 0) then
            za = real(-iza,kind=wp)
            rFinal(:,ipa,ipb,1) = rKappa(:)*(Txyz(:,1,ixa,ixb)*Rnxyz(:,2,iya,iyb)* &
                                             (Two*Rnxyz(:,3,iza+1,izb)*Alpha(:)+za*Rnxyz(:,3,iza-1,izb))+ &
                                             Rnxyz(:,1,ixa,ixb)*Txyz(:,2,iya,iyb)* &
                                             (Two*Rnxyz(:,3,iza+1,izb)*Alpha(:)+za*Rnxyz(:,3,iza-1,izb))+ &
                                             Rnxyz(:,1,ixa,ixb)*Rnxyz(:,2,iya,iyb)* &
                                             (Two*Txyz(:,3,iza+1,izb)*Alpha(:)+za*Txyz(:,3,iza-1,izb)))
          else
            rFinal(:,ipa,ipb,1) = rKappa(:)*Two*Alpha(:)*(Txyz(:,1,ixa,ixb)*Rnxyz(:,2,iya,iyb)*Rnxyz(:,3,iza+1,izb)+ &
                                                          Rnxyz(:,1,ixa,ixb)*Txyz(:,2,iya,iyb)*Rnxyz(:,3,iza+1,izb)+ &
                                                          Rnxyz(:,1,ixa,ixb)*Rnxyz(:,2,iya,iyb)*Txyz(:,3,iza+1,izb))
          end if
        end if
        if (IfGrad(3,2)) then
          if (izb > 0) then
            zb = real(-izb,kind=wp)
            rFinal(:,ipa,ipb,1) = rKappa(:)*(Txyz(:,1,ixa,ixb)*Rnxyz(:,2,iya,iyb)* &
                                             (Two*Rnxyz(:,3,iza,izb+1)*Beta(:)+zb*Rnxyz(:,3,iza,izb-1))+ &
                                             Rnxyz(:,1,ixa,ixb)*Txyz(:,2,iya,iyb)* &
                                             (Two*Rnxyz(:,3,iza,izb+1)*Beta(:)+zb*Rnxyz(:,3,iza,izb-1))+ &
                                             Rnxyz(:,1,ixa,ixb)*Rnxyz(:,2,iya,iyb)* &
                                             (Two*Txyz(:,3,iza,izb+1)*Beta(:)+zb*Txyz(:,3,iza,izb-1)))
          else
            rFinal(:,ipa,ipb,1) = rKappa(:)*Two*Beta(:)*(Txyz(:,1,ixa,ixb)*Rnxyz(:,2,iya,iyb)*Rnxyz(:,3,iza,izb+1)+ &
                                                         Rnxyz(:,1,ixa,ixb)*Txyz(:,2,iya,iyb)*Rnxyz(:,3,iza,izb+1)+ &
                                                         Rnxyz(:,1,ixa,ixb)*Rnxyz(:,2,iya,iyb)*Txyz(:,3,iza,izb+1))
          end if
        end if

      end do
    end do
  end do
end do

return

end subroutine CmbnT1_mck
