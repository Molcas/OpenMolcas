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

use Index_Functions, only: C_Ind, nTri_Elem1
use Constants, only: Two, OneHalf
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: nZeta, la, lb
real(kind=wp), intent(in) :: Rnxyz(nZeta,3,0:la+1,0:lb+1), Zeta(nZeta), Alpha(nZeta), Beta(nZeta)
real(kind=wp), intent(inout) :: rKappa(nZeta)
real(kind=wp), intent(out) :: rFinal(nZeta,nTri_Elem1(la),nTri_Elem1(lb),1)
logical(kind=iwp), intent(in) :: IfGrad(3,2)
integer(kind=iwp) :: ipa, ipb, ixa, ixb, iya, iyaMax, iyb, iybMax, iza, izb
real(kind=wp) :: xa, xb, ya, yb, za, zb

!iRout = 134
!iPrint = nPrint(iRout)

!ii = nTri3_Elem(la)
!jj = nTri3_Elem(lb)
rKappa(:) = rKappa*Zeta**(-OneHalf)
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
            !rFinal(:,ipa,ipb,1) = papb*rKappa(:)* &
            rFinal(:,ipa,ipb,1) = rKappa(:)* &
                                  (Two*Alpha(:)*Rnxyz(:,1,ixa+1,ixb)+xa*Rnxyz(:,1,ixa-1,ixb))*Rnxyz(:,2,iya,iyb)*Rnxyz(:,3,iza,izb)
          else
            !rFinal(:,ipa,ipb,1) = papb*rKappa(:)*Two*Alpha(:)*Rnxyz(:,1,ixa+1,ixb)*Rnxyz(:,2,iya,iyb)*Rnxyz(:,3,iza,izb)
            rFinal(:,ipa,ipb,1) = rKappa(:)*Two*Alpha(:)*Rnxyz(:,1,ixa+1,ixb)*Rnxyz(:,2,iya,iyb)*Rnxyz(:,3,iza,izb)
          end if
        end if
        if (IfGrad(1,2)) then
          if (ixb > 0) then
            xb = real(-ixb,kind=wp)
            !rFinal(:,ipa,ipb,1) = papb*rKappa(:)* &
            rFinal(:,ipa,ipb,1) = rKappa(:)* &
                                  (Two*Beta(:)*Rnxyz(:,1,ixa,ixb+1)+xb*Rnxyz(:,1,ixa,ixb-1))*Rnxyz(:,2,iya,iyb)*Rnxyz(:,3,iza,izb)
          else
            !rFinal(:,ipa,ipb,1) = papb*rKappa(:)*Two*Beta(:)*Rnxyz(:,1,ixa,ixb+1)*Rnxyz(:,2,iya,iyb)*Rnxyz(:,3,iza,izb)
            rFinal(:,ipa,ipb,1) = rKappa(:)*Two*Beta(:)*Rnxyz(:,1,ixa,ixb+1)*Rnxyz(:,2,iya,iyb)*Rnxyz(:,3,iza,izb)
          end if
        end if
        if (IfGrad(2,1)) then
          if (iya > 0) then
            ya = real(-iya,kind=wp)
            !rFinal(:,ipa,ipb,1) = papb*rKappa(:)* &
            rFinal(:,ipa,ipb,1) = rKappa(:)* &
                                  Rnxyz(:,1,ixa,ixb)*(Two*Alpha(:)*Rnxyz(:,2,iya+1,iyb)+ya*Rnxyz(:,2,iya-1,iyb))*Rnxyz(:,3,iza,izb)
          else
            !rFinal(:,ipa,ipb,1) = papb*rKappa(:)*Rnxyz(:,1,ixa,ixb)*Two*Alpha(:)*Rnxyz(:,2,iya+1,iyb)*Rnxyz(:,3,iza,izb)
            rFinal(:,ipa,ipb,1) = rKappa(:)*Rnxyz(:,1,ixa,ixb)*Two*Alpha(:)*Rnxyz(:,2,iya+1,iyb)*Rnxyz(:,3,iza,izb)
          end if
        end if
        if (IfGrad(2,2)) then
          if (iyb > 0) then
            yb = real(-iyb,kind=wp)
            !rFinal(:,ipa,ipb,1) = papb*rKappa(:)* &
            rFinal(:,ipa,ipb,1) = rKappa(:)* &
                                  Rnxyz(:,1,ixa,ixb)*(Two*Beta(:)*Rnxyz(:,2,iya,iyb+1)+yb*Rnxyz(:,2,iya,iyb-1))*Rnxyz(:,3,iza,izb)
          else
            !rFinal(:,ipa,ipb,1) = papb*rKappa(:)*Rnxyz(:,1,ixa,ixb)*Two*Beta(:)*Rnxyz(:,2,iya,iyb+1)*Rnxyz(:,3,iza,izb)
            rFinal(:,ipa,ipb,1) = rKappa(:)*Rnxyz(:,1,ixa,ixb)*Two*Beta(:)*Rnxyz(:,2,iya,iyb+1)*Rnxyz(:,3,iza,izb)
          end if
        end if
        if (IfGrad(3,1)) then
          if (iza > 0) then
            za = real(-iza,kind=wp)
            !rFinal(:,ipa,ipb,1) = papb*rKappa(:)* &
            rFinal(:,ipa,ipb,1) = rKappa(:)* &
                                  Rnxyz(:,1,ixa,ixb)*Rnxyz(:,2,iya,iyb)*(Two*Alpha(:)*Rnxyz(:,3,iza+1,izb)+za*Rnxyz(:,3,iza-1,izb))
          else
            !rFinal(:,ipa,ipb,1) = papb*rKappa(:)*Rnxyz(:,1,ixa,ixb)*Rnxyz(:,2,iya,iyb)*Two*Alpha(:)*Rnxyz(:,3,iza+1,izb)
            rFinal(:,ipa,ipb,1) = rKappa(:)*Rnxyz(:,1,ixa,ixb)*Rnxyz(:,2,iya,iyb)*Two*Alpha(:)*Rnxyz(:,3,iza+1,izb)
          end if
        end if
        if (IfGrad(3,2)) then
          if (izb > 0) then
            zb = real(-izb,kind=wp)
            !rFinal(:,ipa,ipb,1) = papb*rKappa(:)* &
            rFinal(:,ipa,ipb,1) = rKappa(:)* &
                                  Rnxyz(:,1,ixa,ixb)*Rnxyz(:,2,iya,iyb)*(Two*Beta(:)*Rnxyz(:,3,iza,izb+1)+zb*Rnxyz(:,3,iza,izb-1))
          else
            !rFinal(:,ipa,ipb,1) = papb*rKappa(:)**Rnxyz(:,1,ixa,ixb)*Rnxyz(:,2,iya,iyb)*Two*Beta(:)*Rnxyz(:,3,iza,izb+1)
            rFinal(:,ipa,ipb,1) = rKappa(:)*Rnxyz(:,1,ixa,ixb)*Rnxyz(:,2,iya,iyb)*Two*Beta(:)*Rnxyz(:,3,iza,izb+1)
          end if
        end if

      end do
    end do
  end do
end do

return

end subroutine CmbnS1_mck
