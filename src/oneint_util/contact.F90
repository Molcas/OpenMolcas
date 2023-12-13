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

subroutine Contact(Zeta,P,nZeta,A,Axyz,la,RB,Bxyz,lb,CoorO,lOper,iChO,nIC,Array,rFinal,iStabM,nStabM,nComp,rKappa)
!***********************************************************************
!                                                                      *
! Object: to compoute the 1-electron contact term.                     *
!                                                                      *
!     Author: Roland Lindh, Dept. Of Theoretical Chemistry,            *
!             University of Lund, Sweden, February '91                 *
!***********************************************************************

use Index_Functions, only: C_Ind, nTri_Elem1
use Constants, only: Zero, One
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: nZeta, la, lb, nComp, lOper(nComp), iCho(nComp), nIC, nStabM, iStabM(0:nStabM-1)
real(kind=wp), intent(in) :: Zeta(nZeta), P(nZeta,3), A(3), RB(3), CoorO(3), rKappa(nZeta)
real(kind=wp), intent(out) :: Axyz(nZeta,3,0:la), Bxyz(nZeta,3,0:lb), Array(nZeta,nTri_Elem1(la),nTri_Elem1(lb))
real(kind=wp), intent(inout) :: rFinal(nZeta,nTri_Elem1(la),nTri_Elem1(lb),nIC)
#include "print.fh"
integer(kind=iwp) :: ia, ib, iCar, iComp, iDCRT(0:7), ipa, ipb, iPrint, iRout, iStabO(0:7), ixa, ixb, iya, iyb, iza, izb, lDCRT, &
                     llOper, LmbdT, nDCRT, nOp, nStabO
real(kind=wp) :: TC(3)
integer(kind=iwp), external :: NrOpr

!                                                                      *
!***********************************************************************
!                                                                      *
iRout = 170
iPrint = nPrint(iRout)
if (iPrint >= 99) then
  call RecPrt(' In Contact: rKappa',' ',rKappa,nZeta,1)
  call RecPrt(' In Contact: Zeta',' ',Zeta,nZeta,1)
  call RecPrt(' In Contact: P',' ',P,nZeta,3)
end if
!                                                                      *
!***********************************************************************
!                                                                      *
llOper = lOper(1)
do iComp=2,nComp
  llOper = ior(llOper,lOper(iComp))
end do
call SOS(iStabO,nStabO,llOper)
call DCR(LmbdT,iStabM,nStabM,iStabO,nStabO,iDCRT,nDCRT)

do lDCRT=0,nDCRT-1
  call OA(iDCRT(lDCRT),CoorO,TC)

  Array(:,:,:) = Zero

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

  ! Modify z-component to carry the the exponential contribution.

  Bxyz(:,3,0) = exp(-Zeta*((TC(1)-P(:,1))**2+(TC(2)-P(:,2))**2+(TC(3)-P(:,3))**2))

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
          Array(:,ipa,ipb) = Array(:,ipa,ipb)+rKappa*Axyz(:,1,ixa)*Axyz(:,2,iya)*Axyz(:,3,iza)*Bxyz(:,1,ixb)*Bxyz(:,2,iyb)* &
                             Bxyz(:,3,izb)
        end do
      end do
    end do
  end do

  ! Accumulate contributions

  nOp = NrOpr(iDCRT(lDCRT))
  call SymAdO(Array,nZeta,la,lb,nComp,rFinal,nIC,nOp,lOper,iChO,One)

end do

return

end subroutine Contact
