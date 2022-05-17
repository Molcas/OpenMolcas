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

subroutine vCff2D(iDum1,iDum2,nRys,Zeta,ZInv,Eta,EInv,nT,Coori,CoorAC,P,Q,la,lb,lc,ld,U2,PAQP,QCPQ,B10,B00,lac,B01)
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

implicit real*8(A-H,O-Z)
#include "real.fh"
#include "print.fh"
real*8 Zeta(nT), ZInv(nT), Eta(nT), EInv(nT), Coori(3,4), CoorAC(3,2), P(nT,3), Q(nT,3), U2(nRys,nT), PAQP(nRys,nT,3), &
       QCPQ(nRys,nT,3), B10(nRys,nT), B00(nRys,nT), B01(nRys,nT)
real*8 tmp
!define _DEBUGPRINT_
#ifdef _DEBUGPRINT_
! Local arrays
character*30 Label
#endif
logical AeqB, CeqD, EQ
#ifdef _DEBUGPRINT_
logical PrintB10, PrintB00, PrintB01

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

nabMax = la+lb
ncdMax = lc+ld
h12 = Half
if ((nabMax >= 2) .and. (ncdMax >= 2)) then
  do iT=1,nT
    do iRys=1,nRys
      tmp = h12*U2(iRys,iT)
      B00(iRys,iT) = tmp
      B10(iRys,iT) = (h12-tmp*Eta(iT))*ZInv(iT)
      B01(iRys,iT) = (h12-tmp*Zeta(iT))*EInv(iT)
    end do
  end do
# ifdef _DEBUGPRINT_
  PrintB10 = .true.
  PrintB01 = .true.
  PrintB00 = .true.
# endif
else if ((ncdMax == 0) .and. (nabMax >= 2)) then
  do iT=1,nT
    do iRys=1,nRys
      B10(iRys,iT) = (h12-h12*U2(iRys,iT)*Eta(iT))*ZInv(iT)
    end do
  end do
# ifdef _DEBUGPRINT_
  PrintB10 = .true.
# endif
else if ((nabMax == 0) .and. (ncdMax >= 2)) then
  do iT=1,nT
    do iRys=1,nRys
      B01(iRys,iT) = (h12-h12*U2(iRys,iT)*Zeta(iT))*EInv(iT)
    end do
  end do
# ifdef _DEBUGPRINT_
  PrintB01 = .true.
# endif
else if ((ncdMax == 1) .and. (nabMax >= 2)) then
  do iT=1,nT
    do iRys=1,nRys
      tmp = h12*U2(iRys,iT)
      B00(iRys,iT) = tmp
      B10(iRys,iT) = (h12-tmp*Eta(iT))*ZInv(iT)
    end do
  end do
# ifdef _DEBUGPRINT_
  PrintB10 = .true.
  PrintB00 = .true.
# endif
else if ((nabMax == 1) .and. (ncdMax >= 2)) then
  do iT=1,nT
    do iRys=1,nRys
      tmp = h12*U2(iRys,iT)
      B00(iRys,iT) = tmp
      B01(iRys,iT) = (h12-tmp*Zeta(iT))*EInv(iT)
    end do
  end do
# ifdef _DEBUGPRINT_
  PrintB01 = .true.
  PrintB00 = .true.
# endif
else if ((nabMax == 1) .and. (ncdMax == 1)) then
  do iT=1,nT
    do iRys=1,nRys
      B00(iRys,iT) = h12*U2(iRys,iT)
    end do
  end do
# ifdef _DEBUGPRINT_
  PrintB00 = .true.
# endif
end if

if ((nabMax /= 0) .and. (ncdMax /= 0)) then
  if ((.not. AeqB) .and. (.not. CeqD)) then
    do iCar=1,3
      do iT=1,nT
        do iRys=1,nRys
          tmp = U2(iRys,iT)*(Q(iT,iCar)-P(iT,iCar))
          PAQP(iRys,iT,iCar) = P(iT,iCar)-CoorAC(iCar,1)+Eta(iT)*tmp
          QCPQ(iRys,iT,iCar) = Q(iT,iCar)-CoorAC(iCar,2)-Zeta(iT)*tmp
        end do
      end do
    end do
  else if (AeqB .and. (.not. CeqD)) then
    do iCar=1,3
      do iT=1,nT
        do iRys=1,nRys
          tmp = U2(iRys,iT)*(Q(iT,iCar)-P(iT,iCar))
          PAQP(iRys,iT,iCar) = Eta(iT)*tmp
          QCPQ(iRys,iT,iCar) = Q(iT,iCar)-CoorAC(iCar,2)-Zeta(iT)*tmp
        end do
      end do
    end do
  else if ((.not. AeqB) .and. CeqD) then
    do iCar=1,3
      do iT=1,nT
        do iRys=1,nRys
          tmp = U2(iRys,iT)*(Q(iT,iCar)-P(iT,iCar))
          PAQP(iRys,iT,iCar) = P(iT,iCar)-CoorAC(iCar,1)+Eta(iT)*tmp
          QCPQ(iRys,iT,iCar) = -Zeta(iT)*tmp
        end do
      end do
    end do
  else
    do iCar=1,3
      do iT=1,nT
        do iRys=1,nRys
          tmp = U2(iRys,iT)*(Q(iT,iCar)-P(iT,iCar))
          PAQP(iRys,iT,iCar) = Eta(iT)*tmp
          QCPQ(iRys,iT,iCar) = -Zeta(iT)*tmp
        end do
      end do
    end do
  end if
else if (nabMax /= 0) then
  if (.not. AeqB) then
    do iCar=1,3
      do iT=1,nT
        do iRys=1,nRys
          PAQP(iRys,iT,iCar) = P(iT,iCar)-CoorAC(iCar,1)+Eta(iT)*(U2(iRys,iT)*(Q(iT,iCar)-P(iT,iCar)))
        end do
      end do
    end do
  else
    do iCar=1,3
      do iT=1,nT
        do iRys=1,nRys
          PAQP(iRys,iT,iCar) = Eta(iT)*(U2(iRys,iT)*(Q(iT,iCar)-P(iT,iCar)))
        end do
      end do
    end do
  end if
else if (ncdMax /= 0) then
  if (.not. CeqD) then
    do iCar=1,3
      do iT=1,nT
        do iRys=1,nRys
          QCPQ(iRys,iT,iCar) = Q(iT,iCar)-CoorAC(iCar,2)-Zeta(iT)*(U2(iRys,iT)*(Q(iT,iCar)-P(iT,iCar)))
        end do
      end do
    end do
  else
    do iCar=1,3
      do iT=1,nT
        do iRys=1,nRys
          QCPQ(iRys,iT,iCar) = -Zeta(iT)*(U2(iRys,iT)*(Q(iT,iCar)-P(iT,iCar)))
        end do
      end do
    end do
  end if
end if
#ifdef _DEBUGPRINT_
if (la+lb > 0) then
  write(Label,'(A)') ' PAQP(x)'
  !call RecPrt(Label,' ',PAQP(1,1,1),nRys,nT)
  write(Label,'(A)') ' PAQP(y)'
  !call RecPrt(Label,' ',PAQP(1,1,2),nRys,nT)
  write(Label,'(A)') ' PAQP(z)'
  !call RecPrt(Label,' ',PAQP(1,1,3),nRys,nT)
end if
if (lc+ld > 0) then
  write(Label,'(A)') ' QCPQ(x)'
  !call RecPrt(Label,' ',QCPQ(1,1,1),nRys,nT)
  write(Label,'(A)') ' QCPQ(y)'
  !call RecPrt(Label,' ',QCPQ(1,1,2),nRys,nT)
  write(Label,'(A)') ' QCPQ(z)'
  !call RecPrt(Label,' ',QCPQ(1,1,3),nRys,nT)
end if
if (PrintB10) then
  write(Label,'(A)') ' B10'
  call RecPrt(Label,' ',B10(1,1),nRys,nT)
end if
if (PrintB00) then
  write(Label,'(A)') ' B00'
  call RecPrt(Label,' ',B00(1,1),nRys,nT)
end if
if (PrintB01) then
  write(Label,'(A)') ' B01'
  call RecPrt(Label,' ',B01(1,1),nRys,nT)
end if
#endif

return
! Avoid unused argument warnings
if (.false.) then
  call Unused_integer(iDum1)
  call Unused_integer(iDum2)
  call Unused_integer(lac)
end if

end subroutine vCff2D
