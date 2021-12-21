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

subroutine LOOP22(KM,ISTOP,IT1,IT2)

use Definitions, only: wp, iwp

implicit none
integer(kind=iwp) :: KM, ISTOP, IT1, IT2
#include "real_guga.fh"
#include "integ.fh"
integer(kind=iwp) :: IDIF, IWAYKM, KM1
real(kind=wp) :: WM0, WP0

ISTOP = 0
KM1 = KM+1
IDIF = IA(J1(KM1))-IA(J2(KM1))
if ((IDIF < -1) .or. (IDIF > 1)) GO TO 55
if (IDIF < 0) then
  GO TO 51
else if (IDIF == 0) then
  GO TO 52
else
  GO TO 53
end if
51 if (IWAY(KM) == 2) GO TO 55
IWAY(KM) = 2
! GE+
if ((K1(IT1+J1(KM1)) == 0) .or. (K2(IT2+J2(KM1)) == 0)) GO TO 55
if (K0F(JM(KM1)) == 0) GO TO 55
J2(KM) = K2(IT2+J2(KM1))
J1(KM) = J2(KM)
ICOUP1(KM) = ICOUP1(KM1)+IY(IT1+J1(KM1),1)
ICOUP(KM) = ICOUP(KM1)+IY(IT2+J2(KM1),2)
COUP(KM) = COUP(KM1)
GO TO 40
52 IWAYKM = IWAY(KM)
GO TO(59,61,62,55),IWAYKM
59 IWAY(KM) = 2
! (HH+,FF+)
if ((K3(IT1+J1(KM1)) == 0) .or. (K3(IT2+J2(KM1)) == 0)) GO TO 61
WM0 = D0
WP0 = D0
if (K2F(JM1(KM1)) == 0) GO TO 161
J2(KM) = K3(IT2+J2(KM1))
J1(KM) = J2(KM)
ICOUP1(KM) = ICOUP1(KM1)+IY(IT1+J1(KM1),3)
ICOUP(KM) = ICOUP(KM1)+IY(IT2+J2(KM1),3)
WM0 = BS4(IB(J2(KM1)))**2
if (K1F(JM(KM1)) == 0) GO TO 162
GO TO 163
161 if (K1F(JM(KM1)) == 0) GO TO 61
J2(KM) = K3(IT2+J2(KM1))
J1(KM) = J2(KM)
ICOUP1(KM) = ICOUP1(KM1)+IY(IT1+J1(KM1),3)
ICOUP(KM) = ICOUP(KM1)+IY(IT2+J2(KM1),3)
163 WP0 = BS3(IB(J2(KM1))+2)**2
162 COUP(KM) = WM0*COUP1(KM1)+WP0*COUP(KM1)
GO TO 40
! GG+
61 IWAY(KM) = 3
if ((K1(IT1+J1(KM1)) == 0) .or. (K1(IT2+J2(KM1)) == 0)) GO TO 62
if (K0F(JM1(KM1)) == 0) GO TO 62
J2(KM) = K1(IT2+J2(KM1))
J1(KM) = J2(KM)
ICOUP1(KM) = ICOUP1(KM1)+IY(IT1+J1(KM1),1)
ICOUP(KM) = ICOUP(KM1)+IY(IT2+J2(KM1),1)
COUP(KM) = COUP1(KM1)
GO TO 40
! EE+
62 IWAY(KM) = 4
if ((K2(IT1+J1(KM1)) == 0) .or. (K2(IT2+J2(KM1)) == 0)) GO TO 55
if (K0F(JM(KM1)) == 0) GO TO 55
J2(KM) = K2(IT2+J2(KM1))
J1(KM) = J2(KM)
ICOUP1(KM) = ICOUP1(KM1)+IY(IT1+J1(KM1),2)
ICOUP(KM) = ICOUP(KM1)+IY(IT2+J2(KM1),2)
COUP(KM) = COUP(KM1)
GO TO 40
53 if (IWAY(KM) == 2) GO TO 55
IWAY(KM) = 2
! EG+
if ((K2(IT1+J1(KM1)) == 0) .or. (K1(IT2+J2(KM1)) == 0)) GO TO 55
if (K0F(JM1(KM1)) == 0) GO TO 55
J2(KM) = K1(IT2+J2(KM1))
J1(KM) = J2(KM)
ICOUP1(KM) = ICOUP1(KM1)+IY(IT1+J1(KM1),2)
ICOUP(KM) = ICOUP(KM1)+IY(IT2+J2(KM1),1)
COUP(KM) = COUP1(KM1)
GO TO 40
55 ISTOP = 1
40 continue

return

end subroutine LOOP22
