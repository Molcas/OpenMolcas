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

implicit real*8(A-H,O-Z)
dimension IOCR(nIOCR), L0(*), L1(*), L2(*), L3(*)
#include "integ.fh"
dimension IOC(55)

IIN = 0
NREF = 0
IJJ = IV0
KM = 1
J2(KM) = IJJ
11 KM = KM+1
IWAY(KM) = 0
12 KM1 = KM-1
if ((L0(J2(KM1)) == 0) .or. (IWAY(KM) >= 1)) GO TO 14
J2(KM) = L0(J2(KM1))
IWAY(KM) = 1
IOC(KM1) = 0
GO TO 20
14 if ((L1(J2(KM1)) == 0) .or. (IWAY(KM) >= 2)) GO TO 15
J2(KM) = L1(J2(KM1))
IWAY(KM) = 2
IOC(KM1) = 1
GO TO 20
15 if ((L2(J2(KM1)) == 0) .or. (IWAY(KM) >= 3)) GO TO 16
J2(KM) = L2(J2(KM1))
IWAY(KM) = 3
IOC(KM1) = 1
GO TO 20
16 if ((L3(J2(KM1)) == 0) .or. (IWAY(KM) >= 4)) GO TO 17
J2(KM) = L3(J2(KM1))
IWAY(KM) = 4
IOC(KM1) = 2
GO TO 20
17 KM = KM-1
if (KM == 1) GO TO 10
GO TO 12
20 if (KM /= LN+1) GO TO 11
NSJ = 1
do I=1,LN
  if (I > LV) GO TO 113
  if (IOC(I) /= 0) GO TO 12
  GO TO 110
113 if (I > NIORB+LV) GO TO 114
  if (IOC(I) /= 2) GO TO 12
  GO TO 110
114 if (IOC(I) == 1) NSJ = MUL(NSJ,NSM(I))
110 continue
end do
if (NSJ /= LSYM) GO TO 12
NREF = NREF+1
do I=1,LN
  if (I <= NIORB+LV) GO TO 111
  IIN = IIN+1
  if (IIN > nIOCR) then
    write(6,*) 'CIall: IIN > nIOCR'
    write(6,*) 'IIN=',IIN
    write(6,*) 'nIOCR=',nIOCR
    call Abend()
  end if
  IOCR(IIN) = IOC(I)
111 continue
end do
GO TO 12

10 continue

return

end subroutine CIALL
