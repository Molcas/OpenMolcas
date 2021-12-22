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

subroutine LOOP6(KM,ISTOP,IT1,IT2)

use Definitions, only: iwp

implicit none
integer(kind=iwp), intent(in) :: KM, IT1, IT2
integer(kind=iwp), intent(out) :: ISTOP
#include "real_guga.fh"
#include "integ.fh"
integer(kind=iwp) :: IDIF, IWAYKM, KM1

ISTOP = 0
KM1 = KM+1
IDIF = IA(J2(KM1))-IA(J1(KM1))
if ((IDIF < 0) .or. (IDIF > 1)) GO TO 55
if (IDIF == 0) GO TO 45
IWAYKM = IWAY(KM)
GO TO(39,41,42,43,44,55),IWAYKM
! CASE I-L AND Q
39 IWAY(KM) = 2
if ((K0(IT1+J1(KM1)) == 0) .or. (K0(IT2+J2(KM1)) == 0)) GO TO 41
COUP(KM) = COUP(KM1)
J1(KM) = K0(IT1+J1(KM1))
J2(KM) = K0(IT2+J2(KM1))
ICOUP1(KM) = ICOUP1(KM1)
ICOUP(KM) = ICOUP(KM1)
GO TO 40
41 IWAY(KM) = 3
if ((K1(IT1+J1(KM1)) == 0) .or. (K1(IT2+J2(KM1)) == 0)) GO TO 42
COUP(KM) = BL1(IB(J1(KM1))+1)*COUP(KM1)
J1(KM) = K1(IT1+J1(KM1))
J2(KM) = K1(IT2+J2(KM1))
ICOUP1(KM) = ICOUP1(KM1)+IY(IT1+J1(KM1),1)
ICOUP(KM) = ICOUP(KM1)+IY(IT2+J2(KM1),1)
GO TO 40
42 IWAY(KM) = 4
if ((K2(IT1+J1(KM1)) == 0) .or. (K2(IT2+J2(KM1)) == 0)) GO TO 43
COUP(KM) = -COUP(KM1)
J1(KM) = K2(IT1+J1(KM1))
J2(KM) = K2(IT2+J2(KM1))
ICOUP1(KM) = ICOUP1(KM1)+IY(IT1+J1(KM1),2)
ICOUP(KM) = ICOUP(KM1)+IY(IT2+J2(KM1),2)
GO TO 40
43 IWAY(KM) = 5
if ((K3(IT1+J1(KM1)) == 0) .or. (K3(IT2+J2(KM1)) == 0)) GO TO 44
COUP(KM) = -COUP(KM1)
J1(KM) = K3(IT1+J1(KM1))
J2(KM) = K3(IT2+J2(KM1))
ICOUP1(KM) = ICOUP1(KM1)+IY(IT1+J1(KM1),3)
ICOUP(KM) = ICOUP(KM1)+IY(IT2+J2(KM1),3)
GO TO 40
44 IWAY(KM) = 6
if ((K1(IT1+J1(KM1)) == 0) .or. (K2(IT2+J2(KM1)) == 0)) GO TO 55
COUP(KM) = COUP(KM1)/IB(J1(KM1))
J1(KM) = K1(IT1+J1(KM1))
J2(KM) = K2(IT2+J2(KM1))
ICOUP1(KM) = ICOUP1(KM1)+IY(IT1+J1(KM1),1)
ICOUP(KM) = ICOUP(KM1)+IY(IT2+J2(KM1),2)
GO TO 40
! CASE M-P AND R
45 IWAYKM = IWAY(KM)
GO TO(59,61,62,63,64,55),IWAYKM
59 IWAY(KM) = 2
if ((K0(IT1+J1(KM1)) == 0) .or. (K0(IT2+J2(KM1)) == 0)) GO TO 61
COUP(KM) = COUP(KM1)
J1(KM) = K0(IT1+J1(KM1))
J2(KM) = K0(IT2+J2(KM1))
ICOUP1(KM) = ICOUP1(KM1)
ICOUP(KM) = ICOUP(KM1)
GO TO 40
61 IWAY(KM) = 3
if ((K1(IT1+J1(KM1)) == 0) .or. (K1(IT2+J2(KM1)) == 0)) GO TO 62
COUP(KM) = -COUP(KM1)
J1(KM) = K1(IT1+J1(KM1))
J2(KM) = K1(IT2+J2(KM1))
ICOUP1(KM) = ICOUP1(KM1)+IY(IT1+J1(KM1),1)
ICOUP(KM) = ICOUP(KM1)+IY(IT2+J2(KM1),1)
GO TO 40
62 IWAY(KM) = 4
if ((K2(IT1+J1(KM1)) == 0) .or. (K2(IT2+J2(KM1)) == 0)) GO TO 63
COUP(KM) = BL2(IB(J1(KM1))+1)*COUP(KM1)
J1(KM) = K2(IT1+J1(KM1))
J2(KM) = K2(IT2+J2(KM1))
ICOUP1(KM) = ICOUP1(KM1)+IY(IT1+J1(KM1),2)
ICOUP(KM) = ICOUP(KM1)+IY(IT2+J2(KM1),2)
GO TO 40
63 IWAY(KM) = 5
if ((K3(IT1+J1(KM1)) == 0) .or. (K3(IT2+J2(KM1)) == 0)) GO TO 64
COUP(KM) = -COUP(KM1)
J1(KM) = K3(IT1+J1(KM1))
J2(KM) = K3(IT2+J2(KM1))
ICOUP1(KM) = ICOUP1(KM1)+IY(IT1+J1(KM1),3)
ICOUP(KM) = ICOUP(KM1)+IY(IT2+J2(KM1),3)
GO TO 40
64 IWAY(KM) = 6
if ((K2(IT1+J1(KM1)) == 0) .or. (K1(IT2+J2(KM1)) == 0)) GO TO 55
COUP(KM) = -COUP(KM1)/(IB(J1(KM1))+2)
J1(KM) = K2(IT1+J1(KM1))
J2(KM) = K1(IT2+J2(KM1))
ICOUP1(KM) = ICOUP1(KM1)+IY(IT1+J1(KM1),2)
ICOUP(KM) = ICOUP(KM1)+IY(IT2+J2(KM1),1)
GO TO 40
55 ISTOP = 1
40 continue

return

end subroutine LOOP6
