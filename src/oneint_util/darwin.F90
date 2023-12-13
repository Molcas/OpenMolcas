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

subroutine Darwin(Zeta,P,nZeta,A,Axyz,la,RB,Bxyz,lb,rFinal,iStabM,nStabM,nComp,rKappa)
!***********************************************************************
!                                                                      *
! Object: to compoute the 1-electron Darwin contact term.              *
!                                                                      *
!     Author: Roland Lindh, Dept. Of Theoretical Chemistry,            *
!             University of Lund, Sweden, February '91                 *
!***********************************************************************

use Basis_Info, only: dbsc, nCnttp
use Center_Info, only: dc
use Index_Functions, only: C_Ind, nTri_Elem1
use Constants, only: Zero, One, Half, Pi, c_in_au
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: nZeta, la, lb, nStabM, iStabM(0:nStabM-1), nComp
real(kind=wp), intent(in) :: Zeta(nZeta), P(nZeta,3), A(3), RB(3), rKappa(nZeta)
real(kind=wp), intent(out) :: Axyz(nZeta,3,0:la), Bxyz(nZeta,3,0:lb), rFinal(nZeta,nTri_Elem1(la),nTri_Elem1(lb),nComp)
#include "print.fh"
integer(kind=iwp) :: ia, ib, iCar, iDCRT(0:7), ipa, ipb, iPrint, iRout, ixa, ixb, iya, iyb, iza, izb, kCnt, kCnttp, kdc, lDCRT, &
                     LmbdT, nDCRT
real(kind=wp) :: C(3), Fact, Factor, TC(3)

iRout = 170
iPrint = nPrint(iRout)
if (iPrint >= 99) then
  call RecPrt(' In Darwin: rKappa',' ',rKappa,nZeta,1)
  call RecPrt(' In Darwin: Zeta',' ',Zeta,nZeta,1)
  call RecPrt(' In Darwin: P',' ',P,nZeta,3)
end if

rFinal(:,:,:,:) = Zero

kdc = 0
do kCnttp=1,nCnttp
  if (dbsc(kCnttp)%Aux .or. dbsc(kCnttp)%ECP .or. dbsc(kcnttp)%Frag) exit
  do kCnt=1,dbsc(kCnttp)%nCntr
    C(1:3) = dbsc(kCnttp)%Coor(1:3,kCnt)

    call DCR(LmbdT,iStabM,nStabM,dc(kdc+kCnt)%iStab,dc(kdc+kCnt)%nStab,iDCRT,nDCRT)
    Fact = real(nStabM,kind=wp)/real(LmbdT,kind=wp)

    do lDCRT=0,nDCRT-1
      call OA(iDCRT(lDCRT),C,TC)

      ! Compute the value of the angular components associated
      ! to the basis functions centered on the first center.

      Axyz(:,:,0) = One

      if (la /= 0) then
        do iCar=1,3

          Axyz(:,iCar,1) = TC(iCar)-A(iCar)

          do ia=2,la
            Axyz(:,iCar,ia) = Axyz(:,iCar,1)*Axyz(:,iCar,ia-1)
          end do

        end do
      end if

      ! Compute the value of the angular components associated to
      ! the basis functions centered on the second center.

      Bxyz(:,:,0) = One

      ! Modify z-component to carry the charge and the exponential contribution.

      Bxyz(:,3,0) = dbsc(kCnttp)%Charge*exp(-Zeta*((TC(1)-P(:,1))**2+(TC(2)-P(:,2))**2+(TC(3)-P(:,3))**2))

      if (lb /= 0) then
        do iCar=1,3

          Bxyz(:,iCar,1) = TC(iCar)-RB(iCar)

          do ib=2,lb
            Bxyz(:,iCar,ib) = Bxyz(:,iCar,1)*Bxyz(:,iCar,ib-1)
          end do
        end do

        ! Modify z-components with the exponential contribution

        do ib=1,lb
          Bxyz(:,3,ib) = Bxyz(:,3,ib)*Bxyz(:,3,0)
        end do
      end if

      ! Combine contributions from the various angular components.

      do ixa=la,0,-1
        do ixb=lb,0,-1
          do iya=la-ixa,0,-1
            iza = la-ixa-iya
            ipa = C_Ind(la,ixa,iza)
            do iyb=lb-ixb,0,-1
              izb = lb-ixb-iyb
              ipb = C_Ind(lb,ixb,izb)
              rFinal(:,ipa,ipb,1) = rFinal(:,ipa,ipb,1)+Fact*Axyz(:,1,ixa)*Axyz(:,2,iya)*Axyz(:,3,iza)*Bxyz(:,1,ixb)* &
                                    Bxyz(:,2,iyb)*Bxyz(:,3,izb)
            end do
          end do
        end do
      end do

    end do
  end do
  kdc = kdc+dbsc(kCnttp)%nCntr
end do

! Factor from operator (pi/(2*c**2), c=137.036 au)

Factor = Pi*Half/c_in_au**2
do ipa=1,nTri_Elem1(la)
  do ipb=1,nTri_Elem1(lb)
    rFinal(:,ipa,ipb,1) = rKappa*Factor*rFinal(:,ipa,ipb,1)
  end do
end do

return

end subroutine Darwin
