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

subroutine XCff2D(nabMax,ncdMax,nRys,Zeta,ZInv,Eta,EInv,nT,Coori,CoorAC,P,Q,la,lb,lc,ld,U2,PAQP,QCPQ,B10,B00,lac,B01)
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

implicit real*8(A-H,O-Z)
#include "real.fh"
#include "print.fh"
real*8 Zeta(nT), ZInv(nT), Coori(3,4), CoorAC(3,2), P(nT,3), Q(nT,3), U2(nRys,nT), PAQP(nRys,nT,3), QCPQ(nRys,nT,3), &
       B10(nRys,nT,3), B00(nRys,nT,3), B01(nRys,nT,3)
! Local arrays
!define _DEBUGPRINT_
#ifdef _DEBUGPRINT_
character*30 Label
#endif
logical AeqB, CeqD, EQ

#include "macros.fh"
unused_var(nabMax)
unused_var(ncdMax)
unused_var(Eta)
unused_var(EInv)

#ifdef _DEBUGPRINT_
call RecPrt(' In XCff2D: Coori',' ',Coori,3,4)
call RecPrt(' In XCff2D: P',' ',P,nT,3)
call RecPrt(' In XCff2D: Q',' ',Q,nT,3)
#endif
AeqB = EQ(Coori(1,1),Coori(1,2))
CeqD = EQ(Coori(1,3),Coori(1,4))

nabMax = la+lb
ncdMax = ld+lc
h12 = Half

! Compute B10, B00, and B01

if (nabMax > 1) then
  do iT=1,nT
    do iRys=1,nRys
      B10(iRys,iT,1) = (h12-h12*U2(iRys,iT))*ZInv(iT)
    end do
  end do
  call dcopy_(nRys*nT,B10(1,1,1),1,B10(1,1,2),1)
  call dcopy_(nRys*nT,B10(1,1,1),1,B10(1,1,3),1)
end if
if (lac /= 0) then
  call dcopy_(nRys*nT,U2(1,1),1,B00(1,1,1),1)
  call dcopy_(nRys*nT,U2(1,1),1,B00(1,1,2),1)
  call dcopy_(nRys*nT,U2(1,1),1,B00(1,1,3),1)
end if
if (ncdMax > 1) then
  do iT=1,nT
    do iRys=1,nRys
      B01(iRys,iT,1) = Two*Zeta(iT)*U2(iRys,iT)
    end do
  end do
  call dcopy_(nRys*nT,B01(1,1,1),1,B01(1,1,2),1)
  call dcopy_(nRys*nT,B01(1,1,1),1,B01(1,1,3),1)
end if

if ((nabMax /= 0) .and. (ncdMax /= 0)) then
  if ((.not. AeqB) .and. CeqD) then
    do iCar=1,3
      do iT=1,nT
        do iRys=1,nRys
          PAQP(iRys,iT,iCar) = P(iT,iCar)-CoorAC(iCar,1)+(U2(iRys,iT)*(Q(iT,iCar)-P(iT,iCar)))
          QCPQ(iRys,iT,iCar) = -Two*Zeta(iT)*(U2(iRys,iT)*(Q(iT,iCar)-P(iT,iCar)))
        end do
      end do
    end do
  else
    do iCar=1,3
      do iT=1,nT
        do iRys=1,nRys
          PAQP(iRys,iT,iCar) = (U2(iRys,iT)*(Q(iT,iCar)-P(iT,iCar)))
          QCPQ(iRys,iT,iCar) = -Two*Zeta(iT)*(U2(iRys,iT)*(Q(iT,iCar)-P(iT,iCar)))
        end do
      end do
    end do
  end if
else if (nabMax /= 0) then
  if (.not. AeqB) then
    do iCar=1,3
      do iT=1,nT
        do iRys=1,nRys
          PAQP(iRys,iT,iCar) = P(iT,iCar)-CoorAC(iCar,1)+U2(iRys,iT)*(Q(iT,iCar)-P(iT,iCar))
        end do
      end do
    end do
  else
    do iCar=1,3
      do iT=1,nT
        do iRys=1,nRys
          PAQP(iRys,iT,iCar) = U2(iRys,iT)*(Q(iT,iCar)-P(iT,iCar))
        end do
      end do
    end do
  end if
else if (ncdMax /= 0) then
  do iCar=1,3
    do iT=1,nT
      do iRys=1,nRys
        QCPQ(iRys,iT,iCar) = Two*Zeta(iT)*U2(iRys,iT)*(P(iT,iCar)-Q(iT,iCar))
      end do
    end do
  end do
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

end subroutine XCff2D
