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

subroutine Cff2D( &
#                define _CALLING_
#                include "cff2d_interface.fh"
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

use Constants, only: Half
use Definitions, only: wp, iwp

implicit none
#include "cff2d_interface.fh"
integer(kind=iwp) :: iCar, iT
real(kind=wp) :: h12
logical(kind=iwp) :: AeqB, CeqD
logical(kind=iwp), external :: EQ

#ifdef _DEBUGPRINT_
call RecPrt('Cff2D: Coori',' ',Coori,3,4)
call RecPrt('Cff2D: U2',' ',U2,nRys,nT)
#endif
AeqB = EQ(Coori(1,1),Coori(1,2))
CeqD = EQ(Coori(1,3),Coori(1,4))

h12 = Half
if ((nabMax /= 0) .and. (ncdMax /= 0)) then
  B00(:,:,1) = h12*U2(:,:)
  do iT=1,nT
    B10(:,iT,1) = (h12-h12*U2(:,iT)*Eta(iT))*ZInv(iT)
    B01(:,iT,1) = (h12-h12*U2(:,iT)*Zeta(iT))*EInv(iT)
  end do
else if ((ncdMax == 0) .and. (nabMax /= 0) .and. (lac == 0)) then
  do iT=1,nT
    B10(:,iT,1) = (h12-h12*U2(:,iT)*Eta(iT))*ZInv(iT)
  end do
else if ((nabMax == 0) .and. (ncdMax /= 0) .and. (lac == 0)) then
  do iT=1,nT
    B01(:,iT,1) = (h12-h12*U2(:,iT)*Zeta(iT))*EInv(iT)
  end do
else if ((ncdMax == 0) .and. (nabMax /= 0)) then
  B00(:,:,1) = h12*U2(:,:)
  do iT=1,nT
    B10(:,iT,1) = (h12-h12*U2(:,iT)*Eta(iT))*ZInv(iT)
  end do
else if ((nabMax == 0) .and. (ncdMax /= 0)) then
  B00(:,:,1) = h12*U2(:,:)
  do iT=1,nT
    B01(:,iT,1) = (h12-h12*U2(:,iT)*Zeta(iT))*EInv(iT)
  end do
else if ((nabMax == 0) .and. (ncdMax == 0) .and. (lac /= 0)) then
  B00(:,:,1) = h12*U2(:,:)
end if
if (nabMax /= 0) then
  B10(:,:,2) = B10(:,:,1)
  B10(:,:,3) = B10(:,:,1)
end if
if (lac /= 0) then
  B00(:,:,2) = B00(:,:,1)
  B00(:,:,3) = B00(:,:,1)
end if
if (ncdMax /= 0) then
  B01(:,:,2) = B01(:,:,1)
  B01(:,:,3) = B01(:,:,1)
end if

if ((la+lb+nOrdOp /= 0) .and. (lc+ld+nOrdOp /= 0)) then
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
else if (la+lb+nOrdOp /= 0) then
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
else if (lc+ld+nOrdOp /= 0) then
  if (.not. CeqD) then
    do iCar=1,3
      do iT=1,nT
        QCPQ(:,iT,iCar) = Q(iT,iCar)-CoorAC(iCar,2)+Zeta(iT)*U2(:,iT)*(P(iT,iCar)-Q(iT,iCar))
      end do
    end do
  else
    do iCar=1,3
      do iT=1,nT
        QCPQ(:,iT,iCar) = Zeta(iT)*U2(:,iT)*(P(iT,iCar)-Q(iT,iCar))
      end do
    end do
  end if
end if
#ifdef _DEBUGPRINT_
if (la+lb+nOrdOp > 0) then
  call RecPrt('Cff2D: PAQP(x)',' ',PAQP(:,:,1),nRys,nT)
  call RecPrt('Cff2D: PAQP(y)',' ',PAQP(:,:,2),nRys,nT)
  call RecPrt('Cff2D: PAQP(z)',' ',PAQP(:,:,3),nRys,nT)
end if
if (lc+ld+nOrdOp > 0) then
  call RecPrt('Cff2D: QCPQ(x)',' ',QCPQ(:,:,1),nRys,nT)
  call RecPrt('Cff2D: QCPQ(y)',' ',QCPQ(:,:,2),nRys,nT)
  call RecPrt('Cff2D: QCPQ(z)',' ',QCPQ(:,:,3),nRys,nT)
end if
if (nabMax /= 0) then
  call RecPrt('Cff2D: B10(x)',' ',B10(:,:,1),nRys,nT)
  call RecPrt('Cff2D: B10(y)',' ',B10(:,:,2),nRys,nT)
  call RecPrt('Cff2D: B10(z)',' ',B10(:,:,3),nRys,nT)
end if
if (lac /= 0) then
  call RecPrt('Cff2D: B00(x)',' ',B00(:,:,1),nRys,nT)
  call RecPrt('Cff2D: B00(y)',' ',B00(:,:,2),nRys,nT)
  call RecPrt('Cff2D: B00(z)',' ',B00(:,:,3),nRys,nT)
end if
if (ncdMax /= 0) then
  call RecPrt('Cff2D: B01(x)',' ',B01(:,:,1),nRys,nT)
  call RecPrt('Cff2D: B01(y)',' ',B01(:,:,2),nRys,nT)
  call RecPrt('Cff2D: B01(z)',' ',B01(:,:,3),nRys,nT)
end if
#endif

return

end subroutine Cff2D
