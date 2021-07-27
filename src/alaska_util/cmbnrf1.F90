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

subroutine CmbnRF1(Rnxyz,nZeta,la,lb,lr,Zeta,rKappa,Final,nComp,Fact,Temp,Alpha,Beta,Grad,nGrad,DAO,IfGrad,IndGrd,iStab,jStab,kOp, &
                   EF)
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

implicit real*8(A-H,O-Z)
#include "print.fh"
#include "real.fh"
integer IndGrd(3,2), kOp(2)
logical IfGrad(3,2)
real*8 Final(nZeta,(la+1)*(la+2)/2,(lb+1)*(lb+2)/2,nComp,6), Zeta(nZeta), rKappa(nZeta), Fact(nZeta), Temp(nZeta), &
       Rnxyz(nZeta,3,0:la+1,0:lb+1,0:lr), Grad(nGrad), DAO(nZeta,(la+1)*(la+2)/2,(lb+1)*(lb+2)/2), Alpha(nZeta), Beta(nZeta), &
       EF(nComp)
! Statement function for Cartesian index
Ind(ixyz,ix,iz) = (ixyz-ix)*(ixyz-ix+1)/2+iz+1
iOff(ixyz) = ixyz*(ixyz+1)*(ixyz+2)/6
nElem(i) = (i+1)*(i+2)/2

iRout = 134
iPrint = nPrint(iRout)
if (iPrint >= 99) then
  call RecPrt(' In CmbnRF1: EF',' ',EF,nComp,1)
end if

tTwo = Two

do iZeta=1,nZeta
  Fact(iZeta) = rKappa(iZeta)*Zeta(iZeta)**(-Three/Two)
end do

! Loop over angular components of the basis set

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

        ! Combine multipole moment integrals

        if (IfGrad(1,1)) then
          do ix=0,lr
            do iy=0,lr-ix
              if (ixa > 0) then
                xa = -ixa
                do iZeta=1,nZeta
                  Temp(iZeta) = Fact(iZeta)* &
             (tTwo*Alpha(iZeta)*Rnxyz(iZeta,1,ixa+1,ixb,ix)+ &
                             xa*Rnxyz(iZeta,1,ixa-1,ixb,ix))* &
                                Rnxyz(iZeta,2,iya,iyb,iy)
                end do
              else
                do iZeta=1,nZeta
                  Temp(iZeta) = Fact(iZeta)* &
              tTwo*Alpha(iZeta)*Rnxyz(iZeta,1,ixa+1,ixb,ix)* &
                                Rnxyz(iZeta,2,iya,iyb,iy)
                end do
              end if

              do ir=ix+iy,lr
                iz = ir-ix-iy
                iComp = Ind(ir,ix,iz)+iOff(ir)
                do iZeta=1,nZeta
                  Final(iZeta,ipa,ipb,iComp,1) = Temp(iZeta)*Rnxyz(iZeta,3,iza,izb,iz)
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
              (tTwo*Beta(iZeta)*Rnxyz(iZeta,1,ixa,ixb+1,ix)+ &
                             xb*Rnxyz(iZeta,1,ixa,ixb-1,ix))* &
                                Rnxyz(iZeta,2,iya,iyb,iy)
                end do
              else
                do iZeta=1,nZeta
                  Temp(iZeta) = Fact(iZeta)* &
               tTwo*Beta(iZeta)*Rnxyz(iZeta,1,ixa,ixb+1,ix)* &
                                Rnxyz(iZeta,2,iya,iyb,iy)
                end do
              end if

              do ir=ix+iy,lr
                iz = ir-ix-iy
                iComp = Ind(ir,ix,iz)+iOff(ir)
                do iZeta=1,nZeta
                  Final(iZeta,ipa,ipb,iComp,4) = Temp(iZeta)*Rnxyz(iZeta,3,iza,izb,iz)
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
             (tTwo*Alpha(iZeta)*Rnxyz(iZeta,2,iya+1,iyb,iy)+ &
                             ya*Rnxyz(iZeta,2,iya-1,iyb,iy))
                end do
              else
                do iZeta=1,nZeta
                  Temp(iZeta) = Fact(iZeta)* &
                                Rnxyz(iZeta,1,ixa,ixb,ix)* &
              tTwo*Alpha(iZeta)*Rnxyz(iZeta,2,iya+1,iyb,iy)
                end do
              end if

              do ir=ix+iy,lr
                iz = ir-ix-iy
                iComp = Ind(ir,ix,iz)+iOff(ir)
                do iZeta=1,nZeta
                  Final(iZeta,ipa,ipb,iComp,2) = Temp(iZeta)*Rnxyz(iZeta,3,iza,izb,iz)
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
              (tTwo*Beta(iZeta)*Rnxyz(iZeta,2,iya,iyb+1,iy)+ &
                             yb*Rnxyz(iZeta,2,iya,iyb-1,iy))
                end do
              else
                do iZeta=1,nZeta
                  Temp(iZeta) = Fact(iZeta)* &
                                Rnxyz(iZeta,1,ixa,ixb,ix)* &
               tTwo*Beta(iZeta)*Rnxyz(iZeta,2,iya,iyb+1,iy)
                end do
              end if

              do ir=ix+iy,lr
                iz = ir-ix-iy
                iComp = Ind(ir,ix,iz)+iOff(ir)
                do iZeta=1,nZeta
                  Final(iZeta,ipa,ipb,iComp,5) = Temp(iZeta)*Rnxyz(iZeta,3,iza,izb,iz)
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
                iComp = Ind(ir,ix,iz)+iOff(ir)
                if (iza > 0) then
                  za = -iza
                  do iZeta=1,nZeta
                    Final(iZeta,ipa,ipb,iComp,3) = Temp(iZeta)* &
                                (tTwo*Alpha(iZeta)*Rnxyz(iZeta,3,iza+1,izb,iz)+ &
                                                za*Rnxyz(iZeta,3,iza-1,izb,iz))
                  end do
                else
                  do iZeta=1,nZeta
                    Final(iZeta,ipa,ipb,iComp,3) = Temp(iZeta)*tTwo*Alpha(iZeta)*Rnxyz(iZeta,3,iza+1,izb,iz)
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
                iComp = Ind(ir,ix,iz)+iOff(ir)
                if (izb > 0) then
                  zb = -izb
                  do iZeta=1,nZeta
                    Final(iZeta,ipa,ipb,iComp,6) = Temp(iZeta)*&
                                 (tTwo*Beta(iZeta)*Rnxyz(iZeta,3,iza,izb+1,iz)+ &
                                                zb*Rnxyz(iZeta,3,iza,izb-1,iz))
                  end do
                else
                  do iZeta=1,nZeta
                    Final(iZeta,ipa,ipb,iComp,6) = Temp(iZeta)*tTwo*Beta(iZeta)*Rnxyz(iZeta,3,iza,izb+1,iz)
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
  call RecPrt('In CmbnRF1: DAO',' ',DAO,nZeta,nElem(la)*nElem(lb))
  call RecPrt('In CmbnRF1: Final',' ',Final,nZeta*nElem(la)*nElem(lb)*nComp,6)
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
          ps = dble(iPrmt(kOp(1),iChBas(1+iCar)))
          Fct = dble(iStab)/dble(nIrrep)
        else
          i1 = iCar+3
          i2 = iCar
          ps = dble(iPrmt(kOp(2),iChBas(1+iCar)))
          Fct = ps*dble(jStab)/dble(nIrrep)
        end if
        if (IndGrd(iCar,iCn) < 0) then
          ! Gradient via the translational invariance.
          Grad(iGrad) = Grad(iGrad)-Fct*EF(iEF)*DDot_(nDAO,DAO,1,Final(1,1,1,iEF,i2),1)
        else
          Grad(iGrad) = Grad(iGrad)+Fct*EF(iEF)*DDot_(nDAO,DAO,1,Final(1,1,1,iEF,i1),1)
        end if
      end if
    end do ! End loop over cartesian components, iCar
  end do   ! End loop over centers, iCn
end do     ! End loop over EF components, iEF

return

end subroutine CmbnRF1
