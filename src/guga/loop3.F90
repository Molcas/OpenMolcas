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

subroutine LOOP3(KM,ISTOP,IT1,IT2)

use Definitions, only: iwp

implicit none
integer(kind=iwp), intent(in) :: KM, IT1, IT2
integer(kind=iwp), intent(out) :: ISTOP
#include "real_guga.fh"
#include "integ.fh"
integer(kind=iwp) :: IDIF, IWAYKM, KM1

ISTOP = 0
! STOP THE LOOP
KM1 = KM+1
IDIF = IA(J1(KM1))-IA(J2(KM1))
if ((IDIF < 0) .or. (IDIF > 1)) GO TO 52
if (IDIF == 0) GO TO 60
IWAYKM = IWAY(KM)
GO TO(49,51,52),IWAYKM
49 IWAY(KM) = 2
! CASE E-F
if ((K2(IT1+J1(KM1)) == 0) .or. (K0(IT2+J2(KM1)) == 0)) GO TO 51
COUP(KM) = COUP(KM1)
J2(KM) = K2(IT1+J1(KM1))
J1(KM) = J2(KM)
ICOUP1(KM) = ICOUP1(KM1)+IY(IT1+J1(KM1),2)
ICOUP(KM) = ICOUP(KM1)
GO TO 40
51 IWAY(KM) = 3
if ((K3(IT1+J1(KM1)) == 0) .or. (K1(IT2+J2(KM1)) == 0)) GO TO 52
COUP(KM) = COUP(KM1)*BS3(IB(J2(KM1))+1)
J2(KM) = K3(IT1+J1(KM1))
J1(KM) = J2(KM)
ICOUP1(KM) = ICOUP1(KM1)+IY(IT1+J1(KM1),3)
ICOUP(KM) = ICOUP(KM1)+IY(IT2+J2(KM1),1)
GO TO 40
! CASE G-H
60 IWAYKM = IWAY(KM)
GO TO(64,65,52),IWAYKM
64 IWAY(KM) = 2
if ((K1(IT1+J1(KM1)) == 0) .or. (K0(IT2+J2(KM1)) == 0)) GO TO 65
COUP(KM) = COUP(KM1)
J2(KM) = K1(IT1+J1(KM1))
J1(KM) = J2(KM)
ICOUP1(KM) = ICOUP1(KM1)+IY(IT1+J1(KM1),1)
ICOUP(KM) = ICOUP(KM1)
GO TO 40
65 IWAY(KM) = 3
if ((K3(IT1+J1(KM1)) == 0) .or. (K2(IT2+J2(KM1)) == 0)) GO TO 52
COUP(KM) = COUP(KM1)*BS4(IB(J2(KM1))+1)
J2(KM) = K3(IT1+J1(KM1))
J1(KM) = J2(KM)
ICOUP1(KM) = ICOUP1(KM1)+IY(IT1+J1(KM1),3)
ICOUP(KM) = ICOUP(KM1)+IY(IT2+J2(KM1),2)
GO TO 40
52 ISTOP = 1
40 continue

return

end subroutine LOOP3
