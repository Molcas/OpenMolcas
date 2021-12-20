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

subroutine LOOP12(KM,ISTOP,IT1,IT2)

implicit real*8(A-H,O-Z)
#include "real_guga.fh"
#include "integ.fh"

ISTOP = 0
KM1 = KM+1
J2F = IPO(J2(KM1))
IDIF = IA(J1(KM1))-IA(J2(KM1))
if ((IDIF < 0) .or. (IDIF > 1)) GO TO 55
if (IDIF == 1) GO TO 51
IWAYKM = IWAY(KM)
GO TO(39,41,42,55),IWAYKM
39 IWAY(KM) = 2
! GC
if ((K1(IT1+J1(KM1)) == 0) .or. (K1(IT2+J2(KM1)) == 0)) GO TO 41
if (K0F(J2F) == 0) GO TO 41
J1(KM) = K1(IT1+J1(KM1))
J2(KM) = K1(IT2+J2(KM1))
ICOUP1(KM) = ICOUP1(KM1)+IY(IT1+J1(KM1),1)
ICOUP(KM) = ICOUP(KM1)+IY(IT2+J2(KM1),1)
COUP(KM) = COUP(KM1)
GO TO 40
! HD
41 IWAY(KM) = 3
if ((K3(IT1+J1(KM1)) == 0) .or. (K3(IT2+J2(KM1)) == 0)) GO TO 42
if (K2F(J2F) == 0) GO TO 42
J1(KM) = K3(IT1+J1(KM1))
J2(KM) = K3(IT2+J2(KM1))
ICOUP1(KM) = ICOUP1(KM1)+IY(IT1+J1(KM1),3)
ICOUP(KM) = ICOUP(KM1)+IY(IT2+J2(KM1),3)
COUP(KM) = BS4(IB(J2(KM1))+1)*BS2(IB(J2(KM1))+1)*COUP(KM1)
GO TO 40
! GA
42 IWAY(KM) = 4
if ((K1(IT1+J1(KM1)) == 0) .or. (K2(IT2+J2(KM1)) == 0)) GO TO 55
if (K0F(J2F) == 0) GO TO 55
J1(KM) = K1(IT1+J1(KM1))
J2(KM) = K2(IT2+J2(KM1))
ICOUP1(KM) = ICOUP1(KM1)+IY(IT1+J1(KM1),1)
ICOUP(KM) = ICOUP(KM1)+IY(IT2+J2(KM1),2)
COUP(KM) = COUP(KM1)
GO TO 40
51 IWAYKM = IWAY(KM)
GO TO(59,61,62,55),IWAYKM
59 IWAY(KM) = 2
! EC
if ((K2(IT1+J1(KM1)) == 0) .or. (K1(IT2+J2(KM1)) == 0)) GO TO 61
if (K0F(J2F) == 0) GO TO 61
J1(KM) = K2(IT1+J1(KM1))
J2(KM) = K1(IT2+J2(KM1))
ICOUP1(KM) = ICOUP1(KM1)+IY(IT1+J1(KM1),2)
ICOUP(KM) = ICOUP(KM1)+IY(IT2+J2(KM1),1)
COUP(KM) = COUP(KM1)
GO TO 40
! EA
61 IWAY(KM) = 3
if ((K2(IT1+J1(KM1)) == 0) .or. (K2(IT2+J2(KM1)) == 0)) GO TO 62
if (K0F(J2F) == 0) GO TO 62
J1(KM) = K2(IT1+J1(KM1))
J2(KM) = K2(IT2+J2(KM1))
ICOUP1(KM) = ICOUP1(KM1)+IY(IT1+J1(KM1),2)
ICOUP(KM) = ICOUP(KM1)+IY(IT2+J2(KM1),2)
COUP(KM) = COUP(KM1)
GO TO 40
! FB
62 IWAY(KM) = 4
if ((K3(IT1+J1(KM1)) == 0) .or. (K3(IT2+J2(KM1)) == 0)) GO TO 55
if (K1F(J2F) == 0) GO TO 55
J1(KM) = K3(IT1+J1(KM1))
J2(KM) = K3(IT2+J2(KM1))
ICOUP1(KM) = ICOUP1(KM1)+IY(IT1+J1(KM1),3)
ICOUP(KM) = ICOUP(KM1)+IY(IT2+J2(KM1),3)
COUP(KM) = BS3(IB(J2(KM1))+1)*BS1(IB(J2(KM1))+1)*COUP(KM1)
GO TO 40
55 ISTOP = 1
40 continue

return

end subroutine LOOP12
