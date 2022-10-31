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

subroutine CmbnS1(Rnxyz,nZeta,la,lb,Zeta,rKappa,rFinal,Alpha,Beta,Grad,nGrad,DAO,IfGrad,IndGrd,iStab,jStab,kOp)
!***********************************************************************
!                                                                      *
! Object: compute the gradient of the overlap matrix.                  *
!                                                                      *
!     Author: Roland Lindh, Dept. of Theoretical Chemistry,            *
!             University of Lund, SWEDEN                               *
!             October '91.                                             *
!***********************************************************************

use Symmetry_Info, only: iChBas, nIrrep
use Index_Functions, only: C_Ind
use Constants, only: Two, Three
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: nZeta, la, lb, nGrad, IndGrd(3,2), iStab, jStab, kOp(2)
real(kind=wp), intent(in) :: Rnxyz(nZeta,3,0:la+1,0:lb+1), Zeta(nZeta), Alpha(nZeta), Beta(nZeta), &
                             DAO(nZeta,(la+1)*(la+2)/2,(lb+1)*(lb+2)/2)
real(kind=wp), intent(inout) :: rKappa(nZeta), Grad(nGrad)
real(kind=wp), intent(out) :: rFinal(nZeta,(la+1)*(la+2)/2,(lb+1)*(lb+2)/2,6)
logical(kind=iwp), intent(in) :: IfGrad(3,2)
integer(kind=iwp) :: i1, i2, iCar, iCn, iGrad, ipa, ipb, iPrint, iRout, ixa, ixb, iya, iyaMax, iyb, iybMax, iza, izb, iZeta, nDAO
real(kind=wp) :: Fact, ps, xa, xb, ya, yb, za, zb
real(kind=wp), parameter :: exp32 = -Three/Two
integer(kind=iwp), external :: iPrmt
real(kind=wp), external :: DDot_
#include "print.fh"
#include "nac.fh"

iRout = 134
iPrint = nPrint(iRout)

!ii = la*(la+1)*(la+2)/6
!jj = lb*(lb+1)*(lb+2)/6
do iZeta=1,nZeta
  rKappa(iZeta) = rKappa(iZeta)*Zeta(iZeta)**exp32
end do
if (iPrint >= 99) then
  call RecPrt(' In CmbnS1: Zeta  ',' ',Zeta,1,nZeta)
  call RecPrt(' In CmbnS1: rKappa',' ',rKappa,1,nZeta)
  call RecPrt(' In CmbnS1: Alpha ',' ',Alpha,1,nZeta)
  call RecPrt(' In CmbnS1: Beta  ',' ',Beta,1,nZeta)
end if
do ixa=0,la
  iyaMax = la-ixa
  do ixb=0,lb
    iybMax = lb-ixb
    do iya=0,iyaMax
      iza = la-ixa-iya
      ipa = C_Ind(la,ixa,iza)
      !iChBs = iChBas(ii+ipa)
      !pa = real(iPrmt(kOp(1),iChBs),kind=wp)
      do iyb=0,iybMax
        izb = lb-ixb-iyb
        ipb = C_Ind(lb,ixb,izb)
        !jChBs = iChBas(jj+ipb)
        !papb = real(iPrmt(kOp(2),jChBs),kind=wp)*pa

        ! Combine overlap integrals

        !write(u6,*) ' papb=', papb
        if (IfGrad(1,1)) then
          if (ixa > 0) then
            xa = real(-ixa,kind=wp)
            do iZeta=1,nZeta
              !rFinal(iZeta,ipa,ipb,1) = papb*rKappa(iZeta)* &
              rFinal(iZeta,ipa,ipb,1) = rKappa(iZeta)* &
                      (Two*Alpha(iZeta)*Rnxyz(iZeta,1,ixa+1,ixb)+ &
                                     xa*Rnxyz(iZeta,1,ixa-1,ixb))* &
                                        Rnxyz(iZeta,2,iya,iyb)* &
                                        Rnxyz(iZeta,3,iza,izb)
            end do
          else
            do iZeta=1,nZeta
              !rFinal(iZeta,ipa,ipb,1) = papb*rKappa(iZeta)* &
              rFinal(iZeta,ipa,ipb,1) = rKappa(iZeta)* &
                       Two*Alpha(iZeta)*Rnxyz(iZeta,1,ixa+1,ixb)* &
                                        Rnxyz(iZeta,2,iya,iyb)* &
                                        Rnxyz(iZeta,3,iza,izb)
            end do
          end if
        end if
        if (IfGrad(1,2)) then
          if (ixb > 0) then
            xb = real(-ixb,kind=wp)
            do iZeta=1,nZeta
              !rFinal(iZeta,ipa,ipb,4) = papb*rKappa(iZeta)* &
              rFinal(iZeta,ipa,ipb,4) = rKappa(iZeta)* &
                       (Two*Beta(iZeta)*Rnxyz(iZeta,1,ixa,ixb+1)+ &
                                     xb*Rnxyz(iZeta,1,ixa,ixb-1))* &
                                        Rnxyz(iZeta,2,iya,iyb)* &
                                        Rnxyz(iZeta,3,iza,izb)
            end do
          else
            do iZeta=1,nZeta
              !rFinal(iZeta,ipa,ipb,4) = papb*rKappa(iZeta)* &
              rFinal(iZeta,ipa,ipb,4) = rKappa(iZeta)* &
                        Two*Beta(iZeta)*Rnxyz(iZeta,1,ixa,ixb+1)* &
                                        Rnxyz(iZeta,2,iya,iyb)* &
                                        Rnxyz(iZeta,3,iza,izb)
            end do
          end if
        end if
        if (IfGrad(2,1)) then
          if (iya > 0) then
            ya = real(-iya,kind=wp)
            do iZeta=1,nZeta
              !rFinal(iZeta,ipa,ipb,2) = papb*rKappa(iZeta)* &
              rFinal(iZeta,ipa,ipb,2) = rKappa(iZeta)* &
                                        Rnxyz(iZeta,1,ixa,ixb)* &
                      (Two*Alpha(iZeta)*Rnxyz(iZeta,2,iya+1,iyb)+ &
                                     ya*Rnxyz(iZeta,2,iya-1,iyb))* &
                                        Rnxyz(iZeta,3,iza,izb)
            end do
          else
            do iZeta=1,nZeta
              !rFinal(iZeta,ipa,ipb,2) = papb*rKappa(iZeta)* &
              rFinal(iZeta,ipa,ipb,2) = rKappa(iZeta)* &
                                        Rnxyz(iZeta,1,ixa,ixb)* &
                       Two*Alpha(iZeta)*Rnxyz(iZeta,2,iya+1,iyb)* &
                                        Rnxyz(iZeta,3,iza,izb)
            end do
          end if
        end if
        if (IfGrad(2,2)) then
          if (iyb > 0) then
            yb = real(-iyb,kind=wp)
            do iZeta=1,nZeta
              !rFinal(iZeta,ipa,ipb,5) = papb*rKappa(iZeta)* &
              rFinal(iZeta,ipa,ipb,5) = rKappa(iZeta)* &
                                        Rnxyz(iZeta,1,ixa,ixb)* &
                       (Two*Beta(iZeta)*Rnxyz(iZeta,2,iya,iyb+1)+ &
                                     yb*Rnxyz(iZeta,2,iya,iyb-1))* &
                                        Rnxyz(iZeta,3,iza,izb)
            end do
          else
            do iZeta=1,nZeta
              !rFinal(iZeta,ipa,ipb,5) = papb*rKappa(iZeta)* &
              rFinal(iZeta,ipa,ipb,5) = rKappa(iZeta)* &
                                        Rnxyz(iZeta,1,ixa,ixb)* &
                        Two*Beta(iZeta)*Rnxyz(iZeta,2,iya,iyb+1)* &
                                        Rnxyz(iZeta,3,iza,izb)
            end do
          end if
        end if
        if (IfGrad(3,1)) then
          if (iza > 0) then
            za = real(-iza,kind=wp)
            do iZeta=1,nZeta
              !rFinal(iZeta,ipa,ipb,3) = papb*rKappa(iZeta)* &
              rFinal(iZeta,ipa,ipb,3) = rKappa(iZeta)* &
                                        Rnxyz(iZeta,1,ixa,ixb)* &
                                        Rnxyz(iZeta,2,iya,iyb)* &
                      (Two*Alpha(iZeta)*Rnxyz(iZeta,3,iza+1,izb)+ &
                                     za*Rnxyz(iZeta,3,iza-1,izb))
            end do
          else
            do iZeta=1,nZeta
              !rFinal(iZeta,ipa,ipb,3) = papb*rKappa(iZeta)* &
              rFinal(iZeta,ipa,ipb,3) = rKappa(iZeta)* &
                                        Rnxyz(iZeta,1,ixa,ixb)* &
                                        Rnxyz(iZeta,2,iya,iyb)* &
                       Two*Alpha(iZeta)*Rnxyz(iZeta,3,iza+1,izb)
            end do
          end if
        end if
        if (IfGrad(3,2)) then
          if (izb > 0) then
            zb = real(-izb,kind=wp)
            do iZeta=1,nZeta
              !rFinal(iZeta,ipa,ipb,6) = papb*rKappa(iZeta)* &
              rFinal(iZeta,ipa,ipb,6) = rKappa(iZeta)* &
                                        Rnxyz(iZeta,1,ixa,ixb)* &
                                        Rnxyz(iZeta,2,iya,iyb)* &
                       (Two*Beta(iZeta)*Rnxyz(iZeta,3,iza,izb+1)+ &
                                     zb*Rnxyz(iZeta,3,iza,izb-1))
            end do
          else
            do iZeta=1,nZeta
              !rFinal(iZeta,ipa,ipb,6) = papb*rKappa(iZeta)* &
              rFinal(iZeta,ipa,ipb,6) = rKappa(iZeta)* &
                                        Rnxyz(iZeta,1,ixa,ixb)* &
                                        Rnxyz(iZeta,2,iya,iyb)* &
                        Two*Beta(iZeta)*Rnxyz(iZeta,3,iza,izb+1)
            end do
          end if
        end if

      end do
    end do
  end do
end do

! Trace the gradient integrals

nDAO = nZeta*(la+1)*(la+2)/2*(lb+1)*(lb+2)/2
if (iPrint >= 99) then
  call RecPrt(' S(1)',' ',rFinal,nDAO,6)
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
        ps = real(iPrmt(kOp(1),iChBas(1+iCar)),kind=wp)
        Fact = real(iStab,kind=wp)/real(nIrrep,kind=wp)
      else
        i1 = iCar+3
        i2 = iCar
        ps = real(iPrmt(kOp(2),iChBas(1+iCar)),kind=wp)
        Fact = ps*real(jStab,kind=wp)/real(nIrrep,kind=wp)
      end if
      !write(u6,*) ' Fact=',Fact,ps
      if (IndGrd(iCar,iCn) < 0) then
        if (isCSF) then
          ! The CSF NADC component is not translation invariant,
          ! in fact, it is based on an antisymmetric matrix, so a + sign is needed here
          Grad(iGrad) = Grad(iGrad)+Fact*DDot_(nDAO,DAO,1,rFinal(1,1,1,i2),1)
        else
          ! Gradient via translational invariance.
          Grad(iGrad) = Grad(iGrad)-Fact*DDot_(nDAO,DAO,1,rFinal(1,1,1,i2),1)
        end if
      else
        Grad(iGrad) = Grad(iGrad)+Fact*DDot_(nDAO,DAO,1,rFinal(1,1,1,i1),1)
      end if
    end if
  end do
end do

return

end subroutine CmbnS1
