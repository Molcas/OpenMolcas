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
! Copyright (C) 1990, Roland Lindh                                     *
!               1990, IBM                                              *
!***********************************************************************

subroutine Cff2DS(nabMax,ncdMax,nRys,Zeta,ZInv,Eta,EInv,nT,Coori,CoorAC,P,Q,la,lb,lc,ld,U2,PAQP,QCPQ,B10,B00,lac,B01,nOrdOp)
!***********************************************************************
!                                                                      *
! Object: to compute the coefficients in the three terms recurrence    *
!         relation of the 2D-integrals.                                *
!                                                                      *
!     Author: Roland Lindh, IBM Almaden Research Center, San Jose, CA  *
!             March 1990                                               *
!***********************************************************************

use Constants, only: Zero, Half
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: nabMax, ncdMax, nRys, nT, la, lb, lc, ld, lac, nOrdOp
real(kind=wp), intent(in) :: Zeta(nT), ZInv(nT), Eta(nT), EInv(nT), Coori(3,4), CoorAC(3,2), P(nT,3), Q(nT,3), U2(nRys,nT)
real(kind=wp), intent(inout) :: PAQP(nRys,nT,3), QCPQ(nRys,nT,3), B10(nRys,nT,3), B00(nRys,nT,3), B01(nRys,nT,3)
integer(kind=iwp) :: iCar, iT
real(kind=wp) :: h12
logical(kind=iwp) :: AeqB, CeqD, EQ

#include "macros.fh"
unused_var(EInv)
unused_var(Eta)
unused_var(Q)

#ifdef _DEBUGPRINT_
call RecPrt('Cff2Ds: Coori',' ',Coori,3,4)
call RecPrt('Cff2Ds: U2',' ',U2,nRys,nT)
call RecPrt('Cff2Ds: Zeta',' ',Zeta,1,nT)
call RecPrt('Cff2Ds: ZInv',' ',ZInv,1,nT)
#endif
AeqB = EQ(Coori(1,1),Coori(1,2))
CeqD = EQ(Coori(1,3),Coori(1,4))

h12 = Half
if ((nabMax /= 0) .and. (ncdMax /= 0)) then
  B00(:,:,1) = h12*U2(:,:)
  do iT=1,nT
    B10(:,iT,1) = (h12-h12*U2(:,iT)*Zeta(iT))*ZInv(iT)
  end do
  B01(:,:,1) = B10(:,:,1)
else if ((ncdMax == 0) .and. (nabMax /= 0) .and. (lac == 0)) then
  call WarningMessage(2,'Cff2DS: ncdMax == 0 .and. nabMax /= 0 .and. lac == 0')
  write(u6,*) 'ncdMax,nabMax,lac=',ncdMax,nabMax,lac
  call Abend()
else if ((nabMax == 0) .and. (ncdMax /= 0) .and. (lac == 0)) then
  call WarningMessage(2,'Cff2DS: nabMax == 0 .and. ncdMax /= 0 .and. lac == 0')
  write(u6,*) 'ncdMax,nabMax,lac=',ncdMax,nabMax,lac
  call Abend()
else if ((ncdMax == 0) .and. (nabMax /= 0)) then
  call WarningMessage(2,'Cff2DS: ncdMax == 0 .and. nabMax /= 0')
  write(u6,*) 'ncdMax,nabMax,lac=',ncdMax,nabMax,lac
  call Abend()
else if ((nabMax == 0) .and. (ncdMax /= 0)) then
  call WarningMessage(2,'Cff2DS: nabMax == 0 .and. ncdMax /= 0')
  write(u6,*) 'ncdMax,nabMax,lac=',ncdMax,nabMax,lac
  call Abend()
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
        PAQP(:,iT,iCar) = P(iT,iCar)-CoorAC(iCar,1)
      end do
    end do
    QCPQ(:,:,:) = PAQP(:,:,:)
  else if (AeqB .and. (.not. CeqD)) then
    call WarningMessage(2,'Cff2DS: AeqB .and. .not.CeqD')
    write(u6,*) 'AeqB,CeqD=',AeqB,CeqD
    call Abend()
  else if ((.not. AeqB) .and. CeqD) then
    call WarningMessage(2,'Cff2DS: .not.AeqB .and. CeqD')
    write(u6,*) 'AeqB,CeqD=',AeqB,CeqD
    call Abend()
  else
    PAQP(:,:,:) = Zero
    QCPQ(:,:,:) = Zero
  end if
else if (la+lb+nOrdOp /= 0) then
  call WarningMessage(2,'Cff2DS: la+lb /= 0')
  write(u6,*) 'la,lb=',la,lb
  call Abend()
else if (lc+ld+nOrdOp /= 0) then
  call WarningMessage(2,'Cff2DS: lc+ld /= 0')
  write(u6,*) 'lc,ld=',lc,ld
  call Abend()
end if
#ifdef _DEBUGPRINT_
if (la+lb+nOrdOp > 0) then
  call RecPrt('Cff2DS: PAQP(x)',' ',PAQP(:,:,1),nRys,nT)
  call RecPrt('Cff2DS: PAQP(y)',' ',PAQP(:,:,2),nRys,nT)
  call RecPrt('Cff2DS: PAQP(z)',' ',PAQP(:,:,3),nRys,nT)
end if
if (lc+ld+nOrdOp > 0) then
  call RecPrt('Cff2DS: QCPQ(x)',' ',QCPQ(:,:,1),nRys,nT)
  call RecPrt('Cff2DS: QCPQ(y)',' ',QCPQ(:,:,2),nRys,nT)
  call RecPrt('Cff2DS: QCPQ(z)',' ',QCPQ(:,:,3),nRys,nT)
end if
if (nabMax /= 0) then
  call RecPrt('Cff2DS: B10(x)',' ',B10(:,:,1),nRys,nT)
  call RecPrt('Cff2DS: B10(y)',' ',B10(:,:,2),nRys,nT)
  call RecPrt('Cff2DS: B10(z)',' ',B10(:,:,3),nRys,nT)
end if
if (lac /= 0) then
  call RecPrt('Cff2DS: B00(x)',' ',B00(:,:,1),nRys,nT)
  call RecPrt('Cff2DS: B00(y)',' ',B00(:,:,2),nRys,nT)
  call RecPrt('Cff2DS: B00(z)',' ',B00(:,:,3),nRys,nT)
end if
if (ncdMax /= 0) then
  call RecPrt('Cff2DS: B01(x)',' ',B01(:,:,1),nRys,nT)
  call RecPrt('Cff2DS: B01(y)',' ',B01(:,:,2),nRys,nT)
  call RecPrt('Cff2DS: B01(z)',' ',B01(:,:,3),nRys,nT)
end if
#endif

return

end subroutine Cff2DS
