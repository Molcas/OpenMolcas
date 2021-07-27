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
! Copyright (C) 1995, Roland Lindh                                     *
!***********************************************************************

subroutine CmbnW1(Welp0,Welm0,Wel0p,Wel0m,nZeta,la,lb,Zeta,rKappa,Final,Alpha,Beta,Grad,nGrad,DAO,IfGrad,IndGrd,iStab,jStab,kOp)
!***********************************************************************
!                                                                      *
! Object: compute the gradient of the Spherical Well integrals         *
!                                                                      *
!     Author: Roland Lindh, Dept. of Theoretical Chemistry,            *
!             University of Lund, SWEDEN                               *
!             May '95.                                                 *
!***********************************************************************

use Symmetry_Info, only: nIrrep, iChBas

implicit real*8(A-H,O-Z)
#include "print.fh"
#include "real.fh"
real*8 Final(nZeta,(la+1)*(la+2)/2,(lb+1)*(lb+2)/2,6), Zeta(nZeta), rKappa(nZeta), Beta(nZeta), &
       Welp0(nZeta,(la+2)*(la+3)/2,(lb+1)*(lb+2)/2), Welm0(nZeta,la*(la+1)/2,(lb+1)*(lb+2)/2), &
       Wel0p(nZeta,(la+1)*(la+2)/2,(lb+2)*(lb+3)/2), Wel0m(nZeta,(la+1)*(la+2)/2,lb*(lb+1)/2), Alpha(nZeta), Grad(nGrad), &
       DAO(nZeta,(la+1)*(la+2)/2,(lb+1)*(lb+2)/2)
logical IfGrad(3,2)
integer IndGrd(3,2), kOp(2)
! Statement function for Cartesian index
Ind(ix,iy,iz) = ix*0+(iy+iz)*(iy+iz+1)/2+iz+1 ! ix*0 added to avoid compiler warning

iRout = 134
iPrint = nPrint(iRout)

if (iPrint >= 99) then
  call RecPrt(' In CmbnW1: Zeta  ',' ',Zeta,1,nZeta)
  call RecPrt(' In CmbnW1: rKappa',' ',rKappa,1,nZeta)
  call RecPrt(' In CmbnW1: Alpha ',' ',Alpha,1,nZeta)
  call RecPrt(' In CmbnW1: Beta  ',' ',Beta,1,nZeta)
  np0 = (la+2)*(la+3)/2*(lb+1)*(lb+2)/2
  nm0 = la*(la+1)/2*(lb+1)*(lb+2)/2
  n0p = (la+1)*(la+2)/2*(lb+2)*(lb+3)/2
  n0m = (la+1)*(la+2)/2*lb*(lb+1)/2
  call RecPrt(' In CmbnW1: Welp0',' ',Welp0,nZeta,np0)
  if (la >= 1) call RecPrt(' In CmbnW1: Welm0',' ',Welm0,nZeta,nm0)
  call RecPrt(' In CmbnW1: Wel0p',' ',Wel0p,nZeta,n0p)
  if (lb >= 1) call RecPrt(' In CmbnW1: Wel0m',' ',Wel0m,nZeta,n0m)
end if

tTwo = Two
nDAO = nZeta*(la+1)*(la+2)/2*(lb+1)*(lb+2)/2
do ixa=0,la
  iyaMax = la-ixa
  do ixb=0,lb
    iybMax = lb-ixb
    do iya=0,iyaMax
      iza = la-ixa-iya
      ipa = Ind(ixa,iya,iza)
      do iyb=0,iybMax
        izb = lb-ixb-iyb
        ipb = Ind(ixb,iyb,izb)

        ! Combine Spherical Well integrals

        if (IfGrad(1,1)) then
          ipp0 = Ind(ixa+1,iya,iza)
          if (ixa > 0) then
            xa = -ixa
            ipm0 = Ind(ixa-1,iya,iza)
            do iZeta=1,nZeta
              Final(iZeta,ipa,ipb,1) = (tTwo*Alpha(iZeta)*Welp0(iZeta,ipp0,ipb)+ &
                                                       xa*Welm0(iZeta,ipm0,ipb))
            end do
          else
            do iZeta=1,nZeta
              Final(iZeta,ipa,ipb,1) = (tTwo*Alpha(iZeta)*Welp0(iZeta,ipp0,ipb))
            end do
          end if
        end if
        if (IfGrad(1,2)) then
          ip0p = Ind(ixb+1,iyb,izb)
          if (ixb > 0) then
            xb = -ixb
            ip0m = Ind(ixb-1,iyb,izb)
            do iZeta=1,nZeta
              Final(iZeta,ipa,ipb,4) = (tTwo*Beta(iZeta)*Wel0p(iZeta,ipa,ip0p)+ &
                                                      xb*Wel0m(iZeta,ipa,ip0m))
            end do
          else
            do iZeta=1,nZeta
              Final(iZeta,ipa,ipb,4) = (tTwo*Beta(iZeta)*Wel0p(iZeta,ipa,ip0p))
            end do
          end if
        end if
        if (IfGrad(2,1)) then
          ipp0 = Ind(ixa,iya+1,iza)
          if (iya > 0) then
            ya = -iya
            ipm0 = Ind(ixa,iya-1,iza)
            do iZeta=1,nZeta
              Final(iZeta,ipa,ipb,2) = (tTwo*Alpha(iZeta)*Welp0(iZeta,ipp0,ipb)+ &
                                                       ya*Welm0(iZeta,ipm0,ipb))
            end do
          else
            do iZeta=1,nZeta
              Final(iZeta,ipa,ipb,2) = (tTwo*Alpha(iZeta)*Welp0(iZeta,ipp0,ipb))
            end do
          end if
        end if
        if (IfGrad(2,2)) then
          ip0p = Ind(ixb,iyb+1,izb)
          if (iyb > 0) then
            yb = -iyb
            ip0m = Ind(ixb,iyb-1,izb)
            do iZeta=1,nZeta
              Final(iZeta,ipa,ipb,5) = (tTwo*Beta(iZeta)*Wel0p(iZeta,ipa,ip0p)+ &
                                                      yb*Wel0m(iZeta,ipa,ip0m))
            end do
          else
            do iZeta=1,nZeta
              Final(iZeta,ipa,ipb,5) = (tTwo*Beta(iZeta)*Wel0p(iZeta,ipa,ip0p))
            end do
          end if
        end if
        if (IfGrad(3,1)) then
          ipp0 = Ind(ixa,iya,iza+1)
          if (iza > 0) then
            za = -iza
            ipm0 = Ind(ixa,iya,iza-1)
            do iZeta=1,nZeta
              Final(iZeta,ipa,ipb,3) = (tTwo*Alpha(iZeta)*Welp0(iZeta,ipp0,ipb)+ &
                                                       za*Welm0(iZeta,ipm0,ipb))
            end do
          else
            do iZeta=1,nZeta
              Final(iZeta,ipa,ipb,3) = (tTwo*Alpha(iZeta)*Welp0(iZeta,ipp0,ipb))
            end do
          end if
        end if
        if (IfGrad(3,2)) then
          ip0p = Ind(ixb,iyb,izb+1)
          if (izb > 0) then
            zb = -izb
            ip0m = Ind(ixb,iyb,izb-1)
            do iZeta=1,nZeta
              Final(iZeta,ipa,ipb,6) = (tTwo*Beta(iZeta)*Wel0p(iZeta,ipa,ip0p)+ &
                                                      zb*Wel0m(iZeta,ipa,ip0m))
            end do
          else
            do iZeta=1,nZeta
              Final(iZeta,ipa,ipb,6) = (tTwo*Beta(iZeta)*Wel0p(iZeta,ipa,ip0p))
            end do
          end if
        end if

      end do
    end do
  end do
end do

! Trace the gradient integrals

if (iPrint >= 99) then
  call RecPrt(' W(1)',' ',Final,nDAO,6)
  call RecPrt('   D ',' ',DAO,nDAO,1)
end if
do iCn=1,2
  do iCar=1,3
    if (IndGrd(iCar,iCn) /= 0) then
      ! Accumulate contribution to the gradient
      iGrad = abs(IndGrd(iCar,iCn))
      if (iCn == 1) then
        i1 = iCar
        i2 = iCar+3
        ps = dble(iPrmt(kOp(1),iChBas(1+iCar)))
        Fact = dble(iStab)/dble(nIrrep)
      else
        i1 = iCar+3
        i2 = iCar
        ps = dble(iPrmt(kOp(2),iChBas(1+iCar)))
        Fact = ps*dble(jStab)/dble(nIrrep)
      end if
      !write(6,*) ' Fact=',Fact,ps
      if (IndGrd(iCar,iCn) < 0) then
        ! Gradient via translational invariance.
        Grad(iGrad) = Grad(iGrad)-Fact*DDot_(nDAO,DAO,1,Final(1,1,1,i2),1)
      else
        Grad(iGrad) = Grad(iGrad)+Fact*DDot_(nDAO,DAO,1,Final(1,1,1,i1),1)
      end if
    end if
  end do
end do

return

end subroutine CmbnW1
