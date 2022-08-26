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

subroutine CmbnEl(Rnxyz,nZeta,la,lb,lr,Zeta,rKappa,rFinal,Fact,Temp,Alpha,Beta,ifgrad,kcar)
!***********************************************************************
!                                                                      *
! Object: to compute gradient integrals for SC Reaction Fields         *
!                                                                      *
!     Author: Roland Lindh, Dept. of Theoretical Chemistry,            *
!             University of Lund, SWEDEN                               *
!             Modified for reaction field calculations July '92        *
!             Modified for gradient calculations May '95               *
!             Modified for trans. prob.   calculations Oct '97         *
!             by Anders Bernhardsson                                   *
!***********************************************************************

use Index_Functions, only: C_Ind, nTri3_Elem, nTri_Elem1
use Constants, only: Two, OneHalf
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp) :: nZeta, la, lb, lr, kcar
real(kind=wp) :: Rnxyz(nZeta,3,0:la+1,0:lb+1,0:lr), Zeta(nZeta), rKappa(nZeta), rFinal(nZeta,nTri_Elem1(la),nTri_Elem1(lb),2), &
                 Fact(nZeta), Temp(nZeta), Alpha(nZeta), Beta(nZeta)
logical(kind=iwp) :: Ifgrad(3,2)
integer(kind=iwp) :: iComp, ipa, ipb, ir, ix, ixa, ixb, iy, iya, iyaMax, iyb, iybMax, iz, iza, izb, iZeta
real(kind=wp) :: xa, xb, ya, yb, za, zb

do iZeta=1,nZeta
  Fact(iZeta) = rKappa(iZeta)*Zeta(iZeta)**(-OneHalf)
end do

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

        if (ifgrad(1,1)) then
          do ix=0,lr
            do iy=0,lr-ix
              if (ixa > 0) then
                xa = -ixa
                do iZeta=1,nZeta
                  Temp(iZeta) = Fact(iZeta)*(Two*Alpha(iZeta)*Rnxyz(iZeta,1,ixa+1,ixb,ix)+xa*Rnxyz(iZeta,1,ixa-1,ixb,ix))* &
                                Rnxyz(iZeta,2,iya,iyb,iy)
                end do
              else
                do iZeta=1,nZeta
                  Temp(iZeta) = Fact(iZeta)*Two*Alpha(iZeta)*Rnxyz(iZeta,1,ixa+1,ixb,ix)*Rnxyz(iZeta,2,iya,iyb,iy)
                end do
              end if

              do ir=ix+iy,lr
                iz = ir-ix-iy
                iComp = C_Ind(ir,ix,iz)+nTri3_Elem(ir)-1
                if (iComp == kcar) then
                  do iZeta=1,nZeta
                    rFinal(iZeta,ipa,ipb,1) = Temp(iZeta)*Rnxyz(iZeta,3,iza,izb,iz)
                  end do
                end if
              end do
            end do
          end do
        end if
        if (ifgrad(1,2)) then
          do ix=0,lr
            do iy=0,lr-ix
              if (ixb > 0) then
                xb = -ixb
                do iZeta=1,nZeta
                  Temp(iZeta) = Fact(iZeta)*(Two*Beta(iZeta)*Rnxyz(iZeta,1,ixa,ixb+1,ix)+xb*Rnxyz(iZeta,1,ixa,ixb-1,ix))* &
                                Rnxyz(iZeta,2,iya,iyb,iy)
                end do
              else
                do iZeta=1,nZeta
                  Temp(iZeta) = Fact(iZeta)*Two*Beta(iZeta)*Rnxyz(iZeta,1,ixa,ixb+1,ix)*Rnxyz(iZeta,2,iya,iyb,iy)
                end do
              end if

              do ir=ix+iy,lr
                iz = ir-ix-iy
                iComp = C_Ind(ir,ix,iz)+nTri3_Elem(ir)-1
                if (iComp == kcar) then
                  do iZeta=1,nZeta
                    rFinal(iZeta,ipa,ipb,2) = Temp(iZeta)*Rnxyz(iZeta,3,iza,izb,iz)
                  end do
                end if
              end do
            end do
          end do
        end if
        if (ifgrad(2,1)) then
          do ix=0,lr
            do iy=0,lr-ix
              if (iya > 0) then
                ya = -iya
                do iZeta=1,nZeta
                  Temp(iZeta) = Fact(iZeta)*Rnxyz(iZeta,1,ixa,ixb,ix)*(Two*Alpha(iZeta)*Rnxyz(iZeta,2,iya+1,iyb,iy)+ &
                                ya*Rnxyz(iZeta,2,iya-1,iyb,iy))
                end do
              else
                do iZeta=1,nZeta
                  Temp(iZeta) = Fact(iZeta)*Rnxyz(iZeta,1,ixa,ixb,ix)*Two*Alpha(iZeta)*Rnxyz(iZeta,2,iya+1,iyb,iy)
                end do
              end if

              do ir=ix+iy,lr
                iz = ir-ix-iy
                iComp = C_Ind(ir,ix,iz)+nTri3_Elem(ir)-1
                if (iComp == kcar) then
                  do iZeta=1,nZeta
                    rFinal(iZeta,ipa,ipb,1) = Temp(iZeta)*Rnxyz(iZeta,3,iza,izb,iz)
                  end do
                end if
              end do
            end do
          end do
        end if
        if (ifgrad(2,2)) then
          do ix=0,lr
            do iy=0,lr-ix
              if (iyb > 0) then
                yb = -iyb
                do iZeta=1,nZeta
                  Temp(iZeta) = Fact(iZeta)*Rnxyz(iZeta,1,ixa,ixb,ix)*(Two*Beta(iZeta)*Rnxyz(iZeta,2,iya,iyb+1,iy)+ &
                                yb*Rnxyz(iZeta,2,iya,iyb-1,iy))
                end do
              else
                do iZeta=1,nZeta
                  Temp(iZeta) = Fact(iZeta)*Rnxyz(iZeta,1,ixa,ixb,ix)*Two*Beta(iZeta)*Rnxyz(iZeta,2,iya,iyb+1,iy)
                end do
              end if

              do ir=ix+iy,lr
                iz = ir-ix-iy
                iComp = C_Ind(ir,ix,iz)+nTri3_Elem(ir)-1
                if (iComp == kcar) then
                  do iZeta=1,nZeta
                    rFinal(iZeta,ipa,ipb,2) = Temp(iZeta)*Rnxyz(iZeta,3,iza,izb,iz)
                  end do
                end if
              end do
            end do
          end do
        end if
        if (ifgrad(3,1)) then
          do ix=0,lr
            do iy=0,lr-ix
              do iZeta=1,nZeta
                Temp(iZeta) = Fact(iZeta)*Rnxyz(iZeta,1,ixa,ixb,ix)*Rnxyz(iZeta,2,iya,iyb,iy)
              end do

              do ir=ix+iy,lr
                iz = ir-ix-iy
                iComp = C_Ind(ir,ix,iz)+nTri3_Elem(ir)-1
                if (iComp == kcar) then
                  if (iza > 0) then
                    za = -iza
                    do iZeta=1,nZeta
                      rFinal(iZeta,ipa,ipb,1) = Temp(iZeta)*(Two*Alpha(iZeta)*Rnxyz(iZeta,3,iza+1,izb,iz)+ &
                                                za*Rnxyz(iZeta,3,iza-1,izb,iz))
                    end do
                  else
                    do iZeta=1,nZeta
                      rFinal(iZeta,ipa,ipb,1) = Temp(iZeta)*Two*Alpha(iZeta)*Rnxyz(iZeta,3,iza+1,izb,iz)
                    end do
                  end if
                end if
              end do
            end do
          end do
        end if
        if (ifgrad(3,2)) then
          do ix=0,lr
            do iy=0,lr-ix
              do iZeta=1,nZeta
                Temp(iZeta) = Fact(iZeta)*Rnxyz(iZeta,1,ixa,ixb,ix)*Rnxyz(iZeta,2,iya,iyb,iy)
              end do

              do ir=ix+iy,lr
                iz = ir-ix-iy
                iComp = C_Ind(ir,ix,iz)+nTri3_Elem(ir)-1
                if (iComp == kcar) then
                  if (izb > 0) then
                    zb = -izb
                    do iZeta=1,nZeta
                      rFinal(iZeta,ipa,ipb,2) = Temp(iZeta)*(Two*Beta(iZeta)*Rnxyz(iZeta,3,iza,izb+1,iz)+ &
                                                zb*Rnxyz(iZeta,3,iza,izb-1,iz))
                    end do
                  else
                    do iZeta=1,nZeta
                      rFinal(iZeta,ipa,ipb,2) = Temp(iZeta)*Two*Beta(iZeta)*Rnxyz(iZeta,3,iza,izb+1,iz)
                    end do
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
