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
! Copyright (C) 1990,1991,1994, Roland Lindh                           *
!               1990, IBM                                              *
!***********************************************************************

subroutine vCff2D(nabMax,ncdMax,nRys,Zeta,ZInv,Eta,EInv,nT,Coori,CoorAC,P,Q,la,lb,lc,ld,U2,PAQP,QCPQ,B10,B00,lac,B01)
!***********************************************************************
!                                                                      *
! Object: to compute the coefficients in the three terms recurrence    *
!         relation of the 2D-integrals.                                *
!                                                                      *
!     Author: Roland Lindh, IBM Almaden Research Center, San Jose, CA  *
!             March '90                                                *
!                                                                      *
!             Modified loop structure for RISC 1991                    *
!             Modified for decreased memory access January '94.        *
!***********************************************************************

use Constants, only: One, Half
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: nabMax, ncdMax, nRys, nT, la, lb, lc, ld, lac
real(kind=wp), intent(in) :: Zeta(nT), ZInv(nT), Eta(nT), EInv(nT), Coori(3,4), CoorAC(3,2), P(nT,3), Q(nT,3), U2(nRys,nT)
real(kind=wp), intent(inout) :: PAQP(nRys,nT,3), QCPQ(nRys,nT,3), B10(nRys,nT), B00(nRys,nT), B01(nRys,nT)
integer(kind=iwp) :: iCar, iT, nabMax_, ncdMax_
logical(kind=iwp) :: AeqB, CeqD
!define _DEBUGPRINT_
#ifdef _DEBUGPRINT_
logical(kind=iwp) :: PrintB00, PrintB01, PrintB10
#endif
logical(kind=iwp), external :: EQ

#include "macros.fh"
unused_var(nabMax)
unused_var(ncdMax)
unused_var(lac)

#ifdef _DEBUGPRINT_
call RecPrt(' In vCff2D: Coori',' ',Coori,3,4)
call RecPrt(' In vCff2D: U2',' ',U2,nRys,nT)
call RecPrt(' in vCff2d: Zeta',' ',Zeta,1,nT)
call RecPrt(' in vCff2d: Eta ',' ',Eta,1,nT)
call RecPrt(' in vCff2d: ZInv',' ',ZInv,1,nT)
call RecPrt(' in vCff2d: EInv',' ',EInv,1,nT)
#endif
AeqB = EQ(Coori(1,1),Coori(1,2))
CeqD = EQ(Coori(1,3),Coori(1,4))
#ifdef _DEBUGPRINT_
PrintB10 = .false.
PrintB01 = .false.
PrintB00 = .false.
#endif

nabMax_ = la+lb
ncdMax_ = lc+ld
if ((nabMax_ >= 2) .and. (ncdMax_ >= 2)) then
  B00(:,:) = Half*U2
  do iT=1,nT
    B10(:,iT) = Half*(One-U2(:,iT)*Eta(iT))*ZInv(iT)
    B01(:,iT) = Half*(One-U2(:,iT)*Zeta(iT))*EInv(iT)
  end do
# ifdef _DEBUGPRINT_
  PrintB10 = .true.
  PrintB01 = .true.
  PrintB00 = .true.
# endif
else if ((ncdMax_ == 0) .and. (nabMax_ >= 2)) then
  do iT=1,nT
    B10(:,iT) = Half*(One-U2(:,iT)*Eta(iT))*ZInv(iT)
  end do
# ifdef _DEBUGPRINT_
  PrintB10 = .true.
# endif
else if ((nabMax_ == 0) .and. (ncdMax_ >= 2)) then
  do iT=1,nT
    B01(:,iT) = Half*(One-U2(:,iT)*Zeta(iT))*EInv(iT)
  end do
# ifdef _DEBUGPRINT_
  PrintB01 = .true.
# endif
else if ((ncdMax_ == 1) .and. (nabMax_ >= 2)) then
  B00(:,:) = Half*U2
  do iT=1,nT
    B10(:,iT) = Half*(One-U2(:,iT)*Eta(iT))*ZInv(iT)
  end do
# ifdef _DEBUGPRINT_
  PrintB10 = .true.
  PrintB00 = .true.
# endif
else if ((nabMax_ == 1) .and. (ncdMax_ >= 2)) then
  B00(:,:) = Half*U2
  do iT=1,nT
    B01(:,iT) = Half*(One-U2(:,iT)*Zeta(iT))*EInv(iT)
  end do
# ifdef _DEBUGPRINT_
  PrintB01 = .true.
  PrintB00 = .true.
# endif
else if ((nabMax_ == 1) .and. (ncdMax_ == 1)) then
  B00(:,:) = Half*U2
# ifdef _DEBUGPRINT_
  PrintB00 = .true.
# endif
end if

if ((nabMax_ /= 0) .and. (ncdMax_ /= 0)) then
  if ((.not. AeqB) .and. (.not. CeqD)) then
    do iCar=1,3
      do iT=1,nT
        PAQP(:,iT,iCar) = P(iT,iCar)-CoorAC(iCar,1)+Eta(iT)*U2(:,iT)*(Q(iT,iCar)-P(iT,iCar))
        QCPQ(:,iT,iCar) = Q(iT,iCar)-CoorAC(iCar,2)-Zeta(iT)*U2(:,iT)*(Q(iT,iCar)-P(iT,iCar))
      end do
    end do
  else if (AeqB .and. (.not. CeqD)) then
    do iCar=1,3
      do iT=1,nT
        PAQP(:,iT,iCar) = Eta(iT)*U2(:,iT)*(Q(iT,iCar)-P(iT,iCar))
        QCPQ(:,iT,iCar) = Q(iT,iCar)-CoorAC(iCar,2)-Zeta(iT)*U2(:,iT)*(Q(iT,iCar)-P(iT,iCar))
      end do
    end do
  else if ((.not. AeqB) .and. CeqD) then
    do iCar=1,3
      do iT=1,nT
        PAQP(:,iT,iCar) = P(iT,iCar)-CoorAC(iCar,1)+Eta(iT)*U2(:,iT)*(Q(iT,iCar)-P(iT,iCar))
        QCPQ(:,iT,iCar) = -Zeta(iT)*U2(:,iT)*(Q(iT,iCar)-P(iT,iCar))
      end do
    end do
  else
    do iCar=1,3
      do iT=1,nT
        PAQP(:,iT,iCar) = Eta(iT)*U2(:,iT)*(Q(iT,iCar)-P(iT,iCar))
        QCPQ(:,iT,iCar) = -Zeta(iT)*U2(:,iT)*(Q(iT,iCar)-P(iT,iCar))
      end do
    end do
  end if
else if (nabMax_ /= 0) then
  if (.not. AeqB) then
    do iCar=1,3
      do iT=1,nT
        PAQP(:,iT,iCar) = P(iT,iCar)-CoorAC(iCar,1)+Eta(iT)*U2(:,iT)*(Q(iT,iCar)-P(iT,iCar))
      end do
    end do
  else
    do iCar=1,3
      do iT=1,nT
        PAQP(:,iT,iCar) = Eta(iT)*U2(:,iT)*(Q(iT,iCar)-P(iT,iCar))
      end do
    end do
  end if
else if (ncdMax_ /= 0) then
  if (.not. CeqD) then
    do iCar=1,3
      do iT=1,nT
        QCPQ(:,iT,iCar) = Q(iT,iCar)-CoorAC(iCar,2)-Zeta(iT)*U2(:,iT)*(Q(iT,iCar)-P(iT,iCar))
      end do
    end do
  else
    do iCar=1,3
      do iT=1,nT
        QCPQ(:,iT,iCar) = -Zeta(iT)*U2(:,iT)*(Q(iT,iCar)-P(iT,iCar))
      end do
    end do
  end if
end if
#ifdef _DEBUGPRINT_
if (la+lb > 0) then
  !call RecPrt(' PAQP(x)',' ',PAQP(:,:,1),nRys,nT)
  !call RecPrt(' PAQP(y)',' ',PAQP(:,:,2),nRys,nT)
  !call RecPrt(' PAQP(z)',' ',PAQP(:,:,3),nRys,nT)
end if
if (lc+ld > 0) then
  !call RecPrt(' QCPQ(x)',' ',QCPQ(:,:,1),nRys,nT)
  !call RecPrt(' QCPQ(y)',' ',QCPQ(:,:,2),nRys,nT)
  !call RecPrt(' QCPQ(z)',' ',QCPQ(:,:,3),nRys,nT)
end if
if (PrintB10) call RecPrt(' B10',' ',B10(:,:),nRys,nT)
if (PrintB00) call RecPrt(' B00',' ',B00(:,:),nRys,nT)
if (PrintB01) call RecPrt(' B01',' ',B01(:,:),nRys,nT)
#endif

return

end subroutine vCff2D
