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

subroutine LOOP26(KM,ISTOP,IFAI,IT1,IT2)

implicit real*8(A-H,O-Z)
#include "real_guga.fh"
#include "integ.fh"

ISTOP = 0
KM1 = KM+1
J2F = IPO(J2(KM1))
! FOUR POSSIBLE STARTS
IWAYKM = IWAY(KM)
GO TO(34,35,36,37,55),IWAYKM
! (CC+,AA+)
34 IWAY(KM) = 2
if ((K0(IT1+J1(KM1)) == 0) .or. (K0(IT2+J2(KM1)) == 0)) GO TO 35
if (K1F(J2F) == 0) GO TO 135
J1(KM) = K0(IT1+J1(KM1))
J2(KM) = K0(IT2+J2(KM1))
ICOUP1(KM) = 0
ICOUP(KM) = 0
COUP1(KM) = D1
JM1(KM) = K1F(J2F)
if (K2F(J2F) == 0) GO TO 40
GO TO 136
135 if (K2F(J2F) == 0) GO TO 137
J1(KM) = K0(IT1+J1(KM1))
J2(KM) = K0(IT2+J2(KM1))
ICOUP1(KM) = 0
ICOUP(KM) = 0
136 COUP(KM) = D1
JM(KM) = K2F(J2F)
GO TO 40
137 if (IFAI == 0) GO TO 35
J1(KM) = K0(IT1+J1(KM1))
J2(KM) = K0(IT2+J2(KM1))
ICOUP1(KM) = 0
ICOUP(KM) = 0
COUP1(KM) = D0
COUP(KM) = D0
GO TO 40
! BB+
35 IWAY(KM) = 3
if ((K1(IT1+J1(KM1)) == 0) .or. (K1(IT2+J2(KM1)) == 0)) GO TO 36
if (K3F(J2F) == 0) GO TO 138
J1(KM) = K1(IT1+J1(KM1))
J2(KM) = K1(IT2+J2(KM1))
JM(KM) = K3F(J2F)
ICOUP1(KM) = IY(IT1+J1(KM1),1)
ICOUP(KM) = IY(IT2+J2(KM1),1)
COUP(KM) = BS1(IB(J2(KM1))+1)**2
GO TO 40
138 if (IFAI == 0) GO TO 36
J1(KM) = K1(IT1+J1(KM1))
J2(KM) = K1(IT2+J2(KM1))
ICOUP1(KM) = IY(IT1+J1(KM1),1)
ICOUP(KM) = IY(IT2+J2(KM1),1)
COUP1(KM) = D0
COUP(KM) = D0
GO TO 40
! DD+
36 IWAY(KM) = 4
if ((K2(IT1+J1(KM1)) == 0) .or. (K2(IT2+J2(KM1)) == 0)) GO TO 37
if (K3F(J2F) == 0) GO TO 139
J1(KM) = K2(IT1+J1(KM1))
J2(KM) = K2(IT2+J2(KM1))
JM1(KM) = K3F(J2F)
ICOUP1(KM) = IY(IT1+J1(KM1),2)
ICOUP(KM) = IY(IT2+J2(KM1),2)
COUP1(KM) = BS2(IB(J2(KM1))+1)**2
GO TO 40
139 if (IFAI == 0) GO TO 37
J1(KM) = K2(IT1+J1(KM1))
J2(KM) = K2(IT2+J2(KM1))
ICOUP1(KM) = IY(IT1+J1(KM1),2)
ICOUP(KM) = IY(IT2+J2(KM1),2)
COUP1(KM) = D0
COUP(KM) = D0
GO TO 40
! BD+
37 IWAY(KM) = 5
if ((K1(IT1+J1(KM1)) == 0) .or. (K2(IT2+J2(KM1)) == 0)) GO TO 55
if (K3F(J2F) == 0) GO TO 55
J1(KM) = K1(IT1+J1(KM1))
J2(KM) = K2(IT2+J2(KM1))
JM1(KM) = K3F(J2F)
ICOUP1(KM) = IY(IT1+J1(KM1),1)
ICOUP(KM) = IY(IT2+J2(KM1),2)
COUP1(KM) = BS1(IB(J2(KM1))+1)*BS2(IB(J2(KM1))+1)
GO TO 40
55 ISTOP = 1
40 continue

return

end subroutine LOOP26
