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

subroutine PATH(KM,ISTOP,IT1,IT2)

use Definitions, only: iwp

implicit none
integer(kind=iwp) :: KM, ISTOP, IT1, IT2
#include "real_guga.fh"
#include "integ.fh"
integer(kind=iwp) :: IWAYKM, KM1

ISTOP = 0
KM1 = KM+1
IWAYKM = IWAY(KM)
GO TO(132,133,134,135,131),IWAYKM
132 IWAY(KM) = 2
if ((K0(IT1+J1(KM1)) == 0) .or. (K0(IT2+J2(KM1)) == 0)) GO TO 133
J2(KM) = K0(IT2+J2(KM1))
J1(KM) = J2(KM)
ICOUP(KM) = ICOUP(KM1)
ICOUP1(KM) = ICOUP1(KM1)
GO TO 140
133 IWAY(KM) = 3
if ((K1(IT1+J1(KM1)) == 0) .or. (K1(IT2+J2(KM1)) == 0)) GO TO 134
J2(KM) = K1(IT2+J2(KM1))
J1(KM) = J2(KM)
ICOUP(KM) = ICOUP(KM1)+IY(IT2+J2(KM1),1)
ICOUP1(KM) = ICOUP1(KM1)+IY(IT1+J1(KM1),1)
GO TO 140
134 IWAY(KM) = 4
if ((K2(IT1+J1(KM1)) == 0) .or. (K2(IT2+J2(KM1)) == 0)) GO TO 135
J2(KM) = K2(IT2+J2(KM1))
J1(KM) = J2(KM)
ICOUP(KM) = ICOUP(KM1)+IY(IT2+J2(KM1),2)
ICOUP1(KM) = ICOUP1(KM1)+IY(IT1+J1(KM1),2)
GO TO 140
135 IWAY(KM) = 5
if ((K3(IT1+J1(KM1)) == 0) .or. (K3(IT2+J2(KM1)) == 0)) GO TO 131
J2(KM) = K3(IT2+J2(KM1))
J1(KM) = J2(KM)
ICOUP(KM) = ICOUP(KM1)+IY(IT2+J2(KM1),3)
ICOUP1(KM) = ICOUP1(KM1)+IY(IT1+J1(KM1),3)
GO TO 140
131 ISTOP = 1
140 continue

return

end subroutine PATH
