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
! Copyright (C) 1990,1991, Roland Lindh                                *
!               1990, IBM                                              *
!***********************************************************************

subroutine XCff2D( &
#                 define _CALLING_
#                 include "cff2d_interface.fh"
                 )
!***********************************************************************
!                                                                      *
! Object: to compute the coefficients in the three terms recurrence    *
!         relation of the 2D-integrals.                                *
!                                                                      *
!     Author: Roland Lindh, IBM Almaden Research Center, San Jose, CA  *
!             March '90                                                *
!                                                                      *
! Modified loop structure for RISC 1991 R. Lindh, Dept. of Theoretical *
! Chemistry, University of Lund, Sweden.                               *
!***********************************************************************

use Constants, only: One, Two, Half
use Definitions, only: wp, iwp

implicit none
#include "cff2d_interface.fh"
integer(kind=iwp) :: iCar, iT, nabMax_, ncdMax_
logical(kind=iwp) :: AeqB, CeqD
logical(kind=iwp), external :: EQ

#include "macros.fh"
unused_var(nabMax)
unused_var(ncdMax)
unused_var(Eta)
unused_var(EInv)

#ifdef _DEBUGPRINT_
call RecPrt('XCff2D: Coori',' ',Coori,3,4)
call RecPrt('XCff2D: P',' ',P,nT,3)
call RecPrt('XCff2D: Q',' ',Q,nT,3)
call RecPrt('XCff2D: Zeta',' ',Zeta,nT,1)
call RecPrt('XCff2D: ZInv',' ',ZInv,nT,1)
call RecPrt('XCff2D: U2',' ',U2,nRys,nT)
#endif
AeqB = EQ(Coori(1,1),Coori(1,2))
CeqD = EQ(Coori(1,3),Coori(1,4))

if (nOrdOp == 0) then
  nabMax_ = la+lb
  ncdMax_ = ld+lc
else
  nabMax_ = la+lb+2
  ncdMax_ = ld+lc+2
end if

! Compute B10, B00, and B01

if (nabMax_ > 1) then
  do iT=1,nT
    B10(:,iT,1) = Half*(One-U2(:,iT))*ZInv(iT)
  end do
  B10(:,:,2) = B10(:,:,1)
  B10(:,:,3) = B10(:,:,1)
end if
if (lac /= 0) then
  B00(:,:,1) = U2(:,:)
  B00(:,:,2) = U2(:,:)
  B00(:,:,3) = U2(:,:)
end if
if (ncdMax_ > 1) then
  do iT=1,nT
    B01(:,iT,1) = Two*Zeta(iT)*U2(:,iT)
  end do
  B01(:,:,2) = B01(:,:,1)
  B01(:,:,3) = B01(:,:,1)
end if

if ((nabMax_ /= 0) .and. (ncdMax_ /= 0)) then
  if ((.not. AeqB) .and. CeqD) then
    do iCar=1,3
      do iT=1,nT
        PAQP(:,iT,iCar) = P(iT,iCar)-CoorAC(iCar,1)+U2(:,iT)*(Q(iT,iCar)-P(iT,iCar))
        QCPQ(:,iT,iCar) = -Two*Zeta(iT)*U2(:,iT)*(Q(iT,iCar)-P(iT,iCar))
      end do
    end do
  else
    do iCar=1,3
      do iT=1,nT
        PAQP(:,iT,iCar) = U2(:,iT)*(Q(iT,iCar)-P(iT,iCar))
        QCPQ(:,iT,iCar) = -Two*Zeta(iT)*U2(:,iT)*(Q(iT,iCar)-P(iT,iCar))
      end do
    end do
  end if
else if (nabMax_ /= 0) then
  if (.not. AeqB) then
    do iCar=1,3
      do iT=1,nT
        PAQP(:,iT,iCar) = P(iT,iCar)-CoorAC(iCar,1)+U2(:,iT)*(Q(iT,iCar)-P(iT,iCar))
      end do
    end do
  else
    do iCar=1,3
      do iT=1,nT
        PAQP(:,iT,iCar) = U2(:,iT)*(Q(iT,iCar)-P(iT,iCar))
      end do
    end do
  end if
else if (ncdMax_ /= 0) then
  do iCar=1,3
    do iT=1,nT
      QCPQ(:,iT,iCar) = Two*Zeta(iT)*U2(:,iT)*(P(iT,iCar)-Q(iT,iCar))
    end do
  end do
end if
#ifdef _DEBUGPRINT_
if (la+lb+nOrdOp > 0) then
  call RecPrt('xCff2D:  PAQP(x)',' ',PAQP(:,:,1),nRys,nT)
  call RecPrt('xCff2D:  PAQP(y)',' ',PAQP(:,:,2),nRys,nT)
  call RecPrt('xCff2D:  PAQP(z)',' ',PAQP(:,:,3),nRys,nT)
end if
if (lc+ld+nOrdOp > 0) then
  call RecPrt('xCff2D:  QCPQ(x)',' ',QCPQ(:,:,1),nRys,nT)
  call RecPrt('xCff2D:  QCPQ(y)',' ',QCPQ(:,:,2),nRys,nT)
  call RecPrt('xCff2D:  QCPQ(z)',' ',QCPQ(:,:,3),nRys,nT)
end if
if (nabMax_ /= 0) then
  call RecPrt('xCff2D:  B10(x)',' ',B10(:,:,1),nRys,nT)
  call RecPrt('xCff2D:  B10(y)',' ',B10(:,:,2),nRys,nT)
  call RecPrt('xCff2D:  B10(z)',' ',B10(:,:,3),nRys,nT)
end if
if (lac /= 0) then
  call RecPrt('xCff2D:  B00(x)',' ',B00(:,:,1),nRys,nT)
  call RecPrt('xCff2D:  B00(y)',' ',B00(:,:,2),nRys,nT)
  call RecPrt('xCff2D:  B00(z)',' ',B00(:,:,3),nRys,nT)
end if
if (ncdMax_ /= 0) then
  call RecPrt('xCff2D:  B01(x)',' ',B01(:,:,1),nRys,nT)
  call RecPrt('xCff2D:  B01(y)',' ',B01(:,:,2),nRys,nT)
  call RecPrt('xCff2D:  B01(z)',' ',B01(:,:,3),nRys,nT)
end if
#endif

return

end subroutine XCff2D
