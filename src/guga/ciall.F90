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
! Copyright (C) 1986, Per E. M. Siegbahn                               *
!***********************************************************************

subroutine CIALL(LSYM,NREF,IOCR,nIOCR,L0,L1,L2,L3,LV)

use guga_global, only: IV0, IWAY, J2, LN, NIORB, NSM
use Symmetry_Info, only: Mul
use Definitions, only: iwp, u6

implicit none
integer(kind=iwp), intent(in) :: LSYM, nIOCR, L0(*), L1(*), L2(*), L3(*), LV
integer(kind=iwp), intent(out) :: NREF
integer(kind=iwp), intent(inout) :: IOCR(nIOCR)
integer(kind=iwp) :: I, IIN, IJJ, IOC(55), KM, KM1, NSJ

IIN = 0
NREF = 0
IJJ = IV0
KM = 1
J2(KM) = IJJ
loop_1: do
  KM = KM+1
  IWAY(KM) = 0
  loop_2: do
    KM1 = KM-1
    if ((L0(J2(KM1)) /= 0) .and. (IWAY(KM) < 1)) then
      J2(KM) = L0(J2(KM1))
      IWAY(KM) = 1
      IOC(KM1) = 0
    else if ((L1(J2(KM1)) /= 0) .and. (IWAY(KM) < 2)) then
      J2(KM) = L1(J2(KM1))
      IWAY(KM) = 2
      IOC(KM1) = 1
    else if ((L2(J2(KM1)) /= 0) .and. (IWAY(KM) < 3)) then
      J2(KM) = L2(J2(KM1))
      IWAY(KM) = 3
      IOC(KM1) = 1
    else if ((L3(J2(KM1)) /= 0) .and. (IWAY(KM) < 4)) then
      J2(KM) = L3(J2(KM1))
      IWAY(KM) = 4
      IOC(KM1) = 2
    else
      KM = KM-1
      if (KM == 1) exit loop_1
      cycle
    end if
    if (KM /= LN+1) exit
    NSJ = 1
    do I=1,LN
      if (I <= LV) then
        if (IOC(I) /= 0) cycle loop_2
      else if (I <= NIORB+LV) then
        if (IOC(I) /= 2) cycle loop_2
      else
        if (IOC(I) == 1) NSJ = Mul(NSJ,NSM(I))
      end if
    end do
    if (NSJ /= LSYM) cycle
    NREF = NREF+1
    do I=1,LN
      if (I <= NIORB+LV) cycle
      IIN = IIN+1
      if (IIN > nIOCR) then
        write(u6,*) 'CIall: IIN > nIOCR'
        write(u6,*) 'IIN=',IIN
        write(u6,*) 'nIOCR=',nIOCR
        call Abend()
      end if
      IOCR(IIN) = IOC(I)
    end do
  end do loop_2
end do loop_1

return

end subroutine CIALL
