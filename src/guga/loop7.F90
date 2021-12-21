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

subroutine LOOP7(KM,ISTOP,IT1,IT2)

use Definitions, only: iwp

implicit none
integer(kind=iwp) :: KM, ISTOP, IT1, IT2
#include "real_guga.fh"
#include "integ.fh"
integer(kind=iwp) :: J2F, KM1

ISTOP = 0
KM1 = KM+1
J2F = IPO(J2(KM1))
if (IWAY(KM) == 2) GO TO 55
IWAY(KM) = 2
! (CB,AD)
if ((K0(IT1+J1(KM1)) == 0) .or. (K3(IT2+J2(KM1)) == 0)) GO TO 55
if (K1F(J2F) == 0) GO TO 141
J1(KM) = K0(IT1+J1(KM1))
J2(KM) = K3(IT2+J2(KM1))
ICOUP1(KM) = 0
ICOUP(KM) = IY(IT2+J2(KM1),3)
COUP1(KM) = BS1(IB(J2(KM1))+1)
JM1(KM) = K1F(J2F)
if (K2F(J2F) == 0) GO TO 40
GO TO 143
141 if (K2F(J2F) == 0) GO TO 55
J1(KM) = K0(IT1+J1(KM1))
J2(KM) = K3(IT2+J2(KM1))
ICOUP1(KM) = 0
ICOUP(KM) = IY(IT2+J2(KM1),3)
143 COUP(KM) = BS2(IB(J2(KM1))+1)
JM(KM) = K2F(J2F)
GO TO 40
55 ISTOP = 1
40 continue

return

end subroutine LOOP7
