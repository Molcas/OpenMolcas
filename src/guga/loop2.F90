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

subroutine LOOP2(KM,ISTOP,IT1,IT2)

use Definitions, only: iwp

implicit none
integer(kind=iwp), intent(in) :: KM, IT1, IT2
integer(kind=iwp), intent(out) :: ISTOP
#include "real_guga.fh"
#include "integ.fh"
integer(kind=iwp) :: IWAYKM, KM1

ISTOP = 0
KM1 = KM+1
! FOUR POSSIBLE STARTS
IWAYKM = IWAY(KM)
GO TO(34,35,36,37,30),IWAYKM
34 IWAY(KM) = 2
if ((K0(IT2+J2(KM1)) == 0) .or. (K2(IT1+J1(KM1)) == 0)) GO TO 35
COUP(KM) = D1
J1(KM) = K2(IT1+J1(KM1))
J2(KM) = K0(IT2+J2(KM1))
ICOUP(KM) = 0
ICOUP1(KM) = IY(IT1+J1(KM1),2)
GO TO 40
35 IWAY(KM) = 3
if ((K1(IT2+J2(KM1)) == 0) .or. (K3(IT1+J1(KM1)) == 0)) GO TO 36
COUP(KM) = BS1(IB(J2(KM1))+1)
J1(KM) = K3(IT1+J1(KM1))
J2(KM) = K1(IT2+J2(KM1))
ICOUP(KM) = IY(IT2+J2(KM1),1)
ICOUP1(KM) = IY(IT1+J1(KM1),3)
GO TO 40
36 IWAY(KM) = 4
if ((K0(IT2+J2(KM1)) == 0) .or. (K1(IT1+J1(KM1)) == 0)) GO TO 37
COUP(KM) = D1
J1(KM) = K1(IT1+J1(KM1))
J2(KM) = K0(IT2+J2(KM1))
ICOUP(KM) = 0
ICOUP1(KM) = IY(IT1+J1(KM1),1)
GO TO 40
37 IWAY(KM) = 5
if ((K2(IT2+J2(KM1)) == 0) .or. (K3(IT1+J1(KM1)) == 0)) GO TO 30
COUP(KM) = BS2(IB(J2(KM1))+1)
J1(KM) = K3(IT1+J1(KM1))
J2(KM) = K2(IT2+J2(KM1))
ICOUP(KM) = IY(IT2+J2(KM1),2)
ICOUP1(KM) = IY(IT1+J1(KM1),3)
GO TO 40
30 ISTOP = 1
40 continue

return

end subroutine LOOP2
