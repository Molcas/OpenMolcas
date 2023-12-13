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
! Copyright (C) 1991,1992,1995, Roland Lindh                           *
!               1997, Anders Bernhardsson                              *
!***********************************************************************

subroutine CmbnEl(Rnxyz,nZeta,la,lb,lr,Zeta,rKappa,rFinal,Fact,Temp,Alpha,Beta,IfGrad,kcar)
!***********************************************************************
!                                                                      *
! Object: to compute gradient integrals for SC Reaction Fields         *
!                                                                      *
!     Author: Roland Lindh, Dept. of Theoretical Chemistry,            *
!             University of Lund, SWEDEN                               *
!             Modified for reaction field calculations July '92        *
!             Modified for gradient calculations May '95               *
!             Modified for trans. prob. calculations Oct '97           *
!             by Anders Bernhardsson                                   *
!***********************************************************************

use Index_Functions, only: C_Ind, nTri3_Elem, nTri_Elem1
use Constants, only: Two, OneHalf
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: nZeta, la, lb, lr, kcar
real(kind=wp), intent(in) :: Rnxyz(nZeta,3,0:la+1,0:lb+1,0:lr), Zeta(nZeta), rKappa(nZeta), Alpha(nZeta), Beta(nZeta)
real(kind=wp), intent(out) :: rFinal(nZeta,nTri_Elem1(la),nTri_Elem1(lb),2), Fact(nZeta), Temp(nZeta)
logical(kind=iwp), intent(in) :: IfGrad(3,2)
integer(kind=iwp) :: iComp, ipa, ipb, ir, ix, ixa, ixb, iy, iya, iyaMax, iyb, iybMax, iz, iza, izb
real(kind=wp) :: xa, xb, ya, yb, za, zb

Fact(:) = rKappa*Zeta**(-OneHalf)

! Loop over angular components of the basis set

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

        ! Combine multipole moment integrals

        if (IfGrad(1,1)) then
          do ix=0,lr
            do iy=0,lr-ix
              if (ixa > 0) then
                xa = -real(ixa,kind=wp)
                Temp(:) = Fact(:)*(Two*Alpha(:)*Rnxyz(:,1,ixa+1,ixb,ix)+xa*Rnxyz(:,1,ixa-1,ixb,ix))*Rnxyz(:,2,iya,iyb,iy)
              else
                Temp(:) = Fact(:)*Two*Alpha(:)*Rnxyz(:,1,ixa+1,ixb,ix)*Rnxyz(:,2,iya,iyb,iy)
              end if

              do ir=ix+iy,lr
                iz = ir-ix-iy
                iComp = C_Ind(ir,ix,iz)+nTri3_Elem(ir)-1
                if (iComp == kcar) rFinal(:,ipa,ipb,1) = Temp(:)*Rnxyz(:,3,iza,izb,iz)
              end do
            end do
          end do
        end if
        if (IfGrad(1,2)) then
          do ix=0,lr
            do iy=0,lr-ix
              if (ixb > 0) then
                xb = -real(ixb,kind=wp)
                Temp(:) = Fact(:)*(Two*Beta(:)*Rnxyz(:,1,ixa,ixb+1,ix)+xb*Rnxyz(:,1,ixa,ixb-1,ix))*Rnxyz(:,2,iya,iyb,iy)
              else
                Temp(:) = Fact(:)*Two*Beta(:)*Rnxyz(:,1,ixa,ixb+1,ix)*Rnxyz(:,2,iya,iyb,iy)
              end if

              do ir=ix+iy,lr
                iz = ir-ix-iy
                iComp = C_Ind(ir,ix,iz)+nTri3_Elem(ir)-1
                if (iComp == kcar) rFinal(:,ipa,ipb,2) = Temp(:)*Rnxyz(:,3,iza,izb,iz)
              end do
            end do
          end do
        end if
        if (IfGrad(2,1)) then
          do ix=0,lr
            do iy=0,lr-ix
              if (iya > 0) then
                ya = -real(iya,kind=wp)
                Temp(:) = Fact(:)*Rnxyz(:,1,ixa,ixb,ix)*(Two*Alpha(:)*Rnxyz(:,2,iya+1,iyb,iy)+ya*Rnxyz(:,2,iya-1,iyb,iy))
              else
                Temp(:) = Fact(:)*Rnxyz(:,1,ixa,ixb,ix)*Two*Alpha(:)*Rnxyz(:,2,iya+1,iyb,iy)
              end if

              do ir=ix+iy,lr
                iz = ir-ix-iy
                iComp = C_Ind(ir,ix,iz)+nTri3_Elem(ir)-1
                if (iComp == kcar) rFinal(:,ipa,ipb,1) = Temp(:)*Rnxyz(:,3,iza,izb,iz)
              end do
            end do
          end do
        end if
        if (IfGrad(2,2)) then
          do ix=0,lr
            do iy=0,lr-ix
              if (iyb > 0) then
                yb = -real(iyb,kind=wp)
                Temp(:) = Fact(:)*Rnxyz(:,1,ixa,ixb,ix)*(Two*Beta(:)*Rnxyz(:,2,iya,iyb+1,iy)+yb*Rnxyz(:,2,iya,iyb-1,iy))
              else
                Temp(:) = Fact(:)*Rnxyz(:,1,ixa,ixb,ix)*Two*Beta(:)*Rnxyz(:,2,iya,iyb+1,iy)
              end if

              do ir=ix+iy,lr
                iz = ir-ix-iy
                iComp = C_Ind(ir,ix,iz)+nTri3_Elem(ir)-1
                if (iComp == kcar) rFinal(:,ipa,ipb,2) = Temp(:)*Rnxyz(:,3,iza,izb,iz)
              end do
            end do
          end do
        end if
        if (IfGrad(3,1)) then
          do ix=0,lr
            do iy=0,lr-ix
              Temp(:) = Fact(:)*Rnxyz(:,1,ixa,ixb,ix)*Rnxyz(:,2,iya,iyb,iy)

              do ir=ix+iy,lr
                iz = ir-ix-iy
                iComp = C_Ind(ir,ix,iz)+nTri3_Elem(ir)-1
                if (iComp == kcar) then
                  if (iza > 0) then
                    za = -real(iza,kind=wp)
                    rFinal(:,ipa,ipb,1) = Temp(:)*(Two*Alpha(:)*Rnxyz(:,3,iza+1,izb,iz)+za*Rnxyz(:,3,iza-1,izb,iz))
                  else
                    rFinal(:,ipa,ipb,1) = Temp(:)*Two*Alpha(:)*Rnxyz(:,3,iza+1,izb,iz)
                  end if
                end if
              end do
            end do
          end do
        end if
        if (IfGrad(3,2)) then
          do ix=0,lr
            do iy=0,lr-ix
              Temp(:) = Fact(:)*Rnxyz(:,1,ixa,ixb,ix)*Rnxyz(:,2,iya,iyb,iy)

              do ir=ix+iy,lr
                iz = ir-ix-iy
                iComp = C_Ind(ir,ix,iz)+nTri3_Elem(ir)-1
                if (iComp == kcar) then
                  if (izb > 0) then
                    zb = -real(izb,kind=wp)
                    rFinal(:,ipa,ipb,2) = Temp(:)*(Two*Beta(:)*Rnxyz(:,3,iza,izb+1,iz)+zb*Rnxyz(:,3,iza,izb-1,iz))
                  else
                    rFinal(:,ipa,ipb,2) = Temp(:)*Two*Beta(:)*Rnxyz(:,3,iza,izb+1,iz)
                  end if
                end if
              end do
            end do
          end do
        end if

      end do
    end do
  end do
end do

return

end subroutine CmbnEl
