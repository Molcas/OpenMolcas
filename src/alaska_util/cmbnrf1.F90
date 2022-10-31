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
!***********************************************************************

subroutine CmbnRF1(Rnxyz,nZeta,la,lb,lr,Zeta,rKappa,rFinal,nComp,Fact,Temp,Alpha,Beta,Grad,nGrad,DAO,IfGrad,IndGrd,iStab,jStab, &
                   kOp,EF)
!***********************************************************************
!                                                                      *
! Object: to compute gradient integrals for SC Reaction Fields         *
!                                                                      *
!     Author: Roland Lindh, Dept. of Theoretical Chemistry,            *
!             University of Lund, SWEDEN                               *
!             Modified for reaction field calculations July '92        *
!             Modified for gradient calculations May '95               *
!***********************************************************************

use Symmetry_Info, only: nIrrep, iChBas
use Index_Functions, only: C_Ind, nTri_Elem1, nTri3_Elem
use Constants, only: Two, Three
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: nZeta, la, lb, lr, nComp, nGrad, IndGrd(3,2), iStab, jStab, kOp(2)
real(kind=wp), intent(in) :: Rnxyz(nZeta,3,0:la+1,0:lb+1,0:lr), Zeta(nZeta), rKappa(nZeta), Alpha(nZeta), Beta(nZeta), &
                             DAO(nZeta,(la+1)*(la+2)/2,(lb+1)*(lb+2)/2), EF(nComp)
real(kind=wp), intent(out) :: rFinal(nZeta,(la+1)*(la+2)/2,(lb+1)*(lb+2)/2,nComp,6), Fact(nZeta), Temp(nZeta)
real(kind=wp), intent(inout) :: Grad(nGrad)
logical(kind=iwp), intent(in) :: IfGrad(3,2)
integer(kind=iwp) :: i1, i2, iCar, iCn, iComp, iEF, iGrad, ipa, ipb, iPrint, ir, iRout, ix, ixa, ixb, iy, iya, iyaMax, iyb, &
                     iybMax, iz, iza, izb, iZeta, nDAO
real(kind=wp) :: Fct, ps, xa, xb, ya, yb, za, zb
real(kind=wp), parameter :: exp32 = -Three/Two
integer(kind=iwp), external :: iPrmt
real(kind=wp), external :: DDot_
#include "print.fh"

iRout = 134
iPrint = nPrint(iRout)
if (iPrint >= 99) then
  call RecPrt(' In CmbnRF1: EF',' ',EF,nComp,1)
end if

do iZeta=1,nZeta
  Fact(iZeta) = rKappa(iZeta)*Zeta(iZeta)**exp32
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

        if (IfGrad(1,1)) then
          do ix=0,lr
            do iy=0,lr-ix
              if (ixa > 0) then
                xa = -ixa
                do iZeta=1,nZeta
                  Temp(iZeta) = Fact(iZeta)* &
              (Two*Alpha(iZeta)*Rnxyz(iZeta,1,ixa+1,ixb,ix)+ &
                             xa*Rnxyz(iZeta,1,ixa-1,ixb,ix))* &
                                Rnxyz(iZeta,2,iya,iyb,iy)
                end do
              else
                do iZeta=1,nZeta
                  Temp(iZeta) = Fact(iZeta)* &
               Two*Alpha(iZeta)*Rnxyz(iZeta,1,ixa+1,ixb,ix)* &
                                Rnxyz(iZeta,2,iya,iyb,iy)
                end do
              end if

              do ir=ix+iy,lr
                iz = ir-ix-iy
                iComp = C_Ind(ir,ix,iz)+nTri3_Elem(ir)
                do iZeta=1,nZeta
                  rFinal(iZeta,ipa,ipb,iComp,1) = Temp(iZeta)*Rnxyz(iZeta,3,iza,izb,iz)
                end do
              end do
            end do
          end do
        end if
        if (IfGrad(1,2)) then
          do ix=0,lr
            do iy=0,lr-ix
              if (ixb > 0) then
                xb = -ixb
                do iZeta=1,nZeta
                  Temp(iZeta) = Fact(iZeta)* &
               (Two*Beta(iZeta)*Rnxyz(iZeta,1,ixa,ixb+1,ix)+ &
                             xb*Rnxyz(iZeta,1,ixa,ixb-1,ix))* &
                                Rnxyz(iZeta,2,iya,iyb,iy)
                end do
              else
                do iZeta=1,nZeta
                  Temp(iZeta) = Fact(iZeta)* &
                Two*Beta(iZeta)*Rnxyz(iZeta,1,ixa,ixb+1,ix)* &
                                Rnxyz(iZeta,2,iya,iyb,iy)
                end do
              end if

              do ir=ix+iy,lr
                iz = ir-ix-iy
                iComp = C_Ind(ir,ix,iz)+nTri3_Elem(ir)
                do iZeta=1,nZeta
                  rFinal(iZeta,ipa,ipb,iComp,4) = Temp(iZeta)*Rnxyz(iZeta,3,iza,izb,iz)
                end do
              end do
            end do
          end do
        end if
        if (IfGrad(2,1)) then
          do ix=0,lr
            do iy=0,lr-ix
              if (iya > 0) then
                ya = -iya
                do iZeta=1,nZeta
                  Temp(iZeta) = Fact(iZeta)* &
                                Rnxyz(iZeta,1,ixa,ixb,ix)* &
              (Two*Alpha(iZeta)*Rnxyz(iZeta,2,iya+1,iyb,iy)+ &
                             ya*Rnxyz(iZeta,2,iya-1,iyb,iy))
                end do
              else
                do iZeta=1,nZeta
                  Temp(iZeta) = Fact(iZeta)* &
                                Rnxyz(iZeta,1,ixa,ixb,ix)* &
               Two*Alpha(iZeta)*Rnxyz(iZeta,2,iya+1,iyb,iy)
                end do
              end if

              do ir=ix+iy,lr
                iz = ir-ix-iy
                iComp = C_Ind(ir,ix,iz)+nTri3_Elem(ir)
                do iZeta=1,nZeta
                  rFinal(iZeta,ipa,ipb,iComp,2) = Temp(iZeta)*Rnxyz(iZeta,3,iza,izb,iz)
                end do
              end do
            end do
          end do
        end if
        if (IfGrad(2,2)) then
          do ix=0,lr
            do iy=0,lr-ix
              if (iyb > 0) then
                yb = -iyb
                do iZeta=1,nZeta
                  Temp(iZeta) = Fact(iZeta)* &
                                Rnxyz(iZeta,1,ixa,ixb,ix)* &
               (Two*Beta(iZeta)*Rnxyz(iZeta,2,iya,iyb+1,iy)+ &
                             yb*Rnxyz(iZeta,2,iya,iyb-1,iy))
                end do
              else
                do iZeta=1,nZeta
                  Temp(iZeta) = Fact(iZeta)* &
                                Rnxyz(iZeta,1,ixa,ixb,ix)* &
                Two*Beta(iZeta)*Rnxyz(iZeta,2,iya,iyb+1,iy)
                end do
              end if

              do ir=ix+iy,lr
                iz = ir-ix-iy
                iComp = C_Ind(ir,ix,iz)+nTri3_Elem(ir)
                do iZeta=1,nZeta
                  rFinal(iZeta,ipa,ipb,iComp,5) = Temp(iZeta)*Rnxyz(iZeta,3,iza,izb,iz)
                end do
              end do
            end do
          end do
        end if
        if (IfGrad(3,1)) then
          do ix=0,lr
            do iy=0,lr-ix
              do iZeta=1,nZeta
                Temp(iZeta) = Fact(iZeta)* &
                              Rnxyz(iZeta,1,ixa,ixb,ix)* &
                              Rnxyz(iZeta,2,iya,iyb,iy)
              end do

              do ir=ix+iy,lr
                iz = ir-ix-iy
                iComp = C_Ind(ir,ix,iz)+nTri3_Elem(ir)
                if (iza > 0) then
                  za = -iza
                  do iZeta=1,nZeta
                    rFinal(iZeta,ipa,ipb,iComp,3) = Temp(iZeta)* &
                                  (Two*Alpha(iZeta)*Rnxyz(iZeta,3,iza+1,izb,iz)+ &
                                                 za*Rnxyz(iZeta,3,iza-1,izb,iz))
                  end do
                else
                  do iZeta=1,nZeta
                    rFinal(iZeta,ipa,ipb,iComp,3) = Temp(iZeta)*Two*Alpha(iZeta)*Rnxyz(iZeta,3,iza+1,izb,iz)
                  end do
                end if
              end do
            end do
          end do
        end if
        if (IfGrad(3,2)) then
          do ix=0,lr
            do iy=0,lr-ix
              do iZeta=1,nZeta
                Temp(iZeta) = Fact(iZeta)* &
                              Rnxyz(iZeta,1,ixa,ixb,ix)* &
                              Rnxyz(iZeta,2,iya,iyb,iy)
              end do

              do ir=ix+iy,lr
                iz = ir-ix-iy
                iComp = C_Ind(ir,ix,iz)+nTri3_Elem(ir)
                if (izb > 0) then
                  zb = -izb
                  do iZeta=1,nZeta
                    rFinal(iZeta,ipa,ipb,iComp,6) = Temp(iZeta)* &
                                   (Two*Beta(iZeta)*Rnxyz(iZeta,3,iza,izb+1,iz)+ &
                                                 zb*Rnxyz(iZeta,3,iza,izb-1,iz))
                  end do
                else
                  do iZeta=1,nZeta
                    rFinal(iZeta,ipa,ipb,iComp,6) = Temp(iZeta)*Two*Beta(iZeta)*Rnxyz(iZeta,3,iza,izb+1,iz)
                  end do
                end if
              end do
            end do
          end do
        end if

      end do
    end do
  end do
end do
if (iPrint >= 99) then
  call RecPrt('In CmbnRF1: DAO',' ',DAO,nZeta,nTri_Elem1(la)*nTri_Elem1(lb))
  call RecPrt('In CmbnRF1: rFinal',' ',rFinal,nZeta*nTri_Elem1(la)*nTri_Elem1(lb)*nComp,6)
end if

! Trace the gradient integrals

nDAO = nZeta*(la+1)*(la+2)/2*(lb+1)*(lb+2)/2
do iEF=1,nComp
  do iCn=1,2
    do iCar=1,3
      if (IndGrd(iCar,iCn) /= 0) then
        ! Accumulate contribution to the gradient
        iGrad = abs(IndGrd(iCar,iCn))
        if (iCn == 1) then
          i1 = iCar
          i2 = iCar+3
          ps = real(iPrmt(kOp(1),iChBas(1+iCar)),kind=wp)
          Fct = real(iStab,kind=wp)/real(nIrrep,kind=wp)
        else
          i1 = iCar+3
          i2 = iCar
          ps = real(iPrmt(kOp(2),iChBas(1+iCar)),kind=wp)
          Fct = ps*real(jStab,kind=wp)/real(nIrrep,kind=wp)
        end if
        if (IndGrd(iCar,iCn) < 0) then
          ! Gradient via the translational invariance.
          Grad(iGrad) = Grad(iGrad)-Fct*EF(iEF)*DDot_(nDAO,DAO,1,rFinal(1,1,1,iEF,i2),1)
        else
          Grad(iGrad) = Grad(iGrad)+Fct*EF(iEF)*DDot_(nDAO,DAO,1,rFinal(1,1,1,iEF,i1),1)
        end if
      end if
    end do ! End loop over cartesian components, iCar
  end do   ! End loop over centers, iCn
end do     ! End loop over EF components, iEF

return

end subroutine CmbnRF1
