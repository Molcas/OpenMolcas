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

subroutine Cff2DS(nabMax,ncdMax,nRys,Zeta,ZInv,Eta,EInv,nT,Coori,CoorAC,P,Q,la,lb,lc,ld,U2,PAQP,QCPQ,B10,B00,lac,B01)
!***********************************************************************
!                                                                      *
! Object: to compute the coefficients in the three terms recurrence    *
!         relation of the 2D-integrals.                                *
!                                                                      *
!     Author: Roland Lindh, IBM Almaden Research Center, San Jose, CA  *
!             March 1990                                               *
!***********************************************************************

implicit real*8(A-H,O-Z)
#include "real.fh"
#include "print.fh"
real*8 Zeta(nT), ZInv(nT), Eta(nT), EInv(nT), Coori(3,4), CoorAC(3,2), P(nT,3), Q(nT,3), U2(nRys,nT), PAQP(nRys,nT,3), &
       QCPQ(nRys,nT,3), B10(nRys,nT,3), B00(nRys,nT,3), B01(nRys,nT,3)
! Local arrays
logical AeqB, CeqD, EQ
!define _DEBUGPRINT_
#ifdef _DEBUGPRINT_
character*30 Label

call RecPrt(' In Cff2Ds: Coori',' ',Coori,3,4)
call RecPrt(' In Cff2Ds: U2',' ',U2,nRys,nT)
#endif
AeqB = EQ(Coori(1,1),Coori(1,2))
CeqD = EQ(Coori(1,3),Coori(1,4))

h12 = Half
if ((nabMax /= 0) .and. (ncdMax /= 0)) then
  do iT=1,nT
    do iRys=1,nRys
      B00(iRys,iT,1) = h12*U2(iRys,iT)
      B10(iRys,iT,1) = (h12-h12*U2(iRys,iT)*Zeta(iT))*ZInv(iT)
      B01(iRys,iT,1) = B10(iRys,iT,1)
    end do
  end do
else if ((ncdMax == 0) .and. (nabMax /= 0) .and. (lac == 0)) then
  call WarningMessage(2,'Cff2DS: ncdMax == 0 .and. nabMax /= 0 .and. lac == 0')
  write(6,*) 'ncdMax,nabMax,lac=',ncdMax,nabMax,lac
  call Abend()
else if ((nabMax == 0) .and. (ncdMax /= 0) .and. (lac == 0)) then
  call WarningMessage(2,'Cff2DS: nabMax == 0 .and. ncdMax /= 0 .and. lac == 0')
  write(6,*) 'ncdMax,nabMax,lac=',ncdMax,nabMax,lac
  call Abend()
else if ((ncdMax == 0) .and. (nabMax /= 0)) then
  call WarningMessage(2,'Cff2DS: ncdMax == 0 .and. nabMax /= 0')
  write(6,*) 'ncdMax,nabMax,lac=',ncdMax,nabMax,lac
  call Abend()
else if ((nabMax == 0) .and. (ncdMax /= 0)) then
  call WarningMessage(2,'Cff2DS: nabMax == 0 .and. ncdMax /= 0')
  write(6,*) 'ncdMax,nabMax,lac=',ncdMax,nabMax,lac
  call Abend()
else if ((nabMax == 0) .and. (ncdMax == 0) .and. (lac /= 0)) then
  call DYaX(nRys*nT,h12,U2(1,1),1,B00(1,1,1),1)
end if
if (nabMax /= 0) then
  call dcopy_(nRys*nT,B10(1,1,1),1,B10(1,1,2),1)
  call dcopy_(nRys*nT,B10(1,1,1),1,B10(1,1,3),1)
end if
if (lac /= 0) then
  call dcopy_(nRys*nT,B00(1,1,1),1,B00(1,1,2),1)
  call dcopy_(nRys*nT,B00(1,1,1),1,B00(1,1,3),1)
end if
if (ncdMax /= 0) then
  call dcopy_(nRys*nT,B01(1,1,1),1,B01(1,1,2),1)
  call dcopy_(nRys*nT,B01(1,1,1),1,B01(1,1,3),1)
end if

if ((la+lb /= 0) .and. (lc+ld /= 0)) then
  if ((.not. AeqB) .and. (.not. CeqD)) then
    do iCar=1,3
      do iT=1,nT
        do iRys=1,nRys
          PAQP(iRys,iT,iCar) = P(iT,iCar)-CoorAC(iCar,1)
          QCPQ(iRys,iT,iCar) = PAQP(iRys,iT,iCar)
        end do
      end do
    end do
  else if (AeqB .and. (.not. CeqD)) then
    call WarningMessage(2,'Cff2DS: AeqB .and. .not.CeqD')
    write(6,*) 'AeqB,CeqD=',AeqB,CeqD
    call Abend()
  else if ((.not. AeqB) .and. CeqD) then
    call WarningMessage(2,'Cff2DS: .not.AeqB .and. CeqD')
    write(6,*) 'AeqB,CeqD=',AeqB,CeqD
    call Abend()
  else
    call dcopy_(3*nRys*nT,[Zero],0,PAQP,1)
    call dcopy_(3*nRys*nT,[Zero],0,QCPQ,1)
  end if
else if (la+lb /= 0) then
  call WarningMessage(2,'Cff2DS: la+lb /= 0')
  write(6,*) 'la,lb=',la,lb
  call Abend()
else if (lc+ld /= 0) then
  call WarningMessage(2,'Cff2DS: lc+ld /= 0')
  write(6,*) 'lc,ld=',lc,ld
  call Abend()
end if
#ifdef _DEBUGPRINT_
if (la+lb > 0) then
  write(Label,'(A)') ' PAQP(x)'
  call RecPrt(Label,' ',PAQP(1,1,1),nRys,nT)
  write(Label,'(A)') ' PAQP(y)'
  call RecPrt(Label,' ',PAQP(1,1,2),nRys,nT)
  write(Label,'(A)') ' PAQP(z)'
  call RecPrt(Label,' ',PAQP(1,1,3),nRys,nT)
end if
if (lc+ld > 0) then
  write(Label,'(A)') ' QCPQ(x)'
  call RecPrt(Label,' ',QCPQ(1,1,1),nRys,nT)
  write(Label,'(A)') ' QCPQ(y)'
  call RecPrt(Label,' ',QCPQ(1,1,2),nRys,nT)
  write(Label,'(A)') ' QCPQ(z)'
  call RecPrt(Label,' ',QCPQ(1,1,3),nRys,nT)
end if
if (nabMax /= 0) then
  write(Label,'(A)') ' B10(x)'
  call RecPrt(Label,' ',B10(1,1,1),nRys,nT)
  write(Label,'(A)') ' B10(y)'
  call RecPrt(Label,' ',B10(1,1,2),nRys,nT)
  write(Label,'(A)') ' B10(z)'
  call RecPrt(Label,' ',B10(1,1,3),nRys,nT)
end if
if (lac /= 0) then
  write(Label,'(A)') ' B00(x)'
  call RecPrt(Label,' ',B00(1,1,1),nRys,nT)
  write(Label,'(A)') ' B00(y)'
  call RecPrt(Label,' ',B00(1,1,2),nRys,nT)
  write(Label,'(A)') ' B00(z)'
  call RecPrt(Label,' ',B00(1,1,3),nRys,nT)
end if
if (ncdMax /= 0) then
  write(Label,'(A)') ' B01(x)'
  call RecPrt(Label,' ',B01(1,1,1),nRys,nT)
  write(Label,'(A)') ' B01(y)'
  call RecPrt(Label,' ',B01(1,1,2),nRys,nT)
  write(Label,'(A)') ' B01(z)'
  call RecPrt(Label,' ',B01(1,1,3),nRys,nT)
end if
#endif

return
! Avoid unused argument warnings
if (.false.) then
  call Unused_real_array(EInv)
  call Unused_real_array(Eta)
  call Unused_real_array(Q)
end if

end subroutine Cff2DS
