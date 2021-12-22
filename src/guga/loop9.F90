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

subroutine LOOP9(KM,ISTOP,IT1,IT2)

use Definitions, only: iwp

implicit none
integer(kind=iwp), intent(in) :: KM, IT1, IT2
integer(kind=iwp), intent(out) :: ISTOP
#include "real_guga.fh"
#include "integ.fh"
integer(kind=iwp) :: IDIF, J1F, KM1

ISTOP = 0
KM1 = KM+1
J1F = IPO(J1(KM1))
IDIF = IA(J1(KM1))-IA(J2(KM1))
if ((IDIF < 0) .or. (IDIF > 1)) GO TO 55
if (IDIF == 1) GO TO 41
if (IWAY(KM) == 2) GO TO 55
IWAY(KM) = 2
! B+G
if ((K3(IT1+J1(KM1)) == 0) .or. (K0(IT2+J2(KM1)) == 0)) GO TO 55
if (K1F(J1F) == 0) GO TO 55
J1(KM) = K3(IT1+J1(KM1))
J2(KM) = K0(IT2+J2(KM1))
ICOUP1(KM) = ICOUP1(KM1)+IY(IT1+J1(KM1),3)
ICOUP(KM) = ICOUP(KM1)
COUP(KM) = BS1(IB(J2(KM1))+2)*COUP(KM1)
GO TO 40
! D+E
41 if (IWAY(KM) == 2) GO TO 55
IWAY(KM) = 2
if ((K3(IT1+J1(KM1)) == 0) .or. (K0(IT2+J2(KM1)) == 0)) GO TO 55
if (K2F(J1F) == 0) GO TO 55
J1(KM) = K3(IT1+J1(KM1))
J2(KM) = K0(IT2+J2(KM1))
ICOUP1(KM) = ICOUP1(KM1)+IY(IT1+J1(KM1),3)
ICOUP(KM) = ICOUP(KM1)
COUP(KM) = BS2(IB(J2(KM1)))*COUP(KM1)
GO TO 40
55 ISTOP = 1
40 continue

return

end subroutine LOOP9
