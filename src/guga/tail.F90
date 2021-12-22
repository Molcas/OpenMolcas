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

subroutine TAIL(LL,IJJ,ITAI,ITAIL,L0,L1,L2,L3,IT1,IT2)

use Definitions, only: iwp

#include "intent.fh"

implicit none
integer(kind=iwp), intent(in) :: LL, IJJ, ITAIL, L0(*), L1(*), L2(*), L3(*), IT1, IT2
integer(kind=iwp), intent(_OUT_) :: ITAI(*)
#include "integ.fh"
integer(kind=iwp) :: I, KM, KM1, L

L = LL+1
if (ITAIL == 0) then
  return
end if
do I=1,ITAIL
  ITAI(I) = 0
end do
if (L == LN+1) ITAI(1) = 1
KM = L
J2(KM) = IJJ
ICOUP(L) = 1
ICOUP1(L) = 1
11 KM = KM+1
IWAY(KM) = 0
12 KM1 = KM-1
if (IWAY(KM) >= 1) GO TO 14
if ((L0(IT1+J2(KM1)) == 0) .or. (L0(IT2+J2(KM1)) == 0)) GO TO 14
J2(KM) = L0(IT2+J2(KM1))
IWAY(KM) = 1
ICOUP(KM) = ICOUP(KM1)
ICOUP1(KM) = ICOUP1(KM1)
GO TO 20
14 if (IWAY(KM) >= 2) GO TO 15
if ((L1(IT1+J2(KM1)) == 0) .or. (L1(IT2+J2(KM1)) == 0)) GO TO 15
J2(KM) = L1(IT2+J2(KM1))
IWAY(KM) = 2
ICOUP(KM) = ICOUP(KM1)+IY(IT2+J2(KM),1)
ICOUP1(KM) = ICOUP1(KM1)+IY(IT1+J2(KM),1)
GO TO 20
15 if (IWAY(KM) >= 3) GO TO 16
if ((L2(IT1+J2(KM1)) == 0) .or. (L2(IT2+J2(KM1)) == 0)) GO TO 16
J2(KM) = L2(IT2+J2(KM1))
IWAY(KM) = 3
ICOUP(KM) = ICOUP(KM1)+IY(IT2+J2(KM),2)
ICOUP1(KM) = ICOUP1(KM1)+IY(IT1+J2(KM),2)
GO TO 20
16 if (IWAY(KM) >= 4) GO TO 17
if ((L3(IT1+J2(KM1)) == 0) .or. (L3(IT2+J2(KM1)) == 0)) GO TO 17
J2(KM) = L3(IT2+J2(KM1))
IWAY(KM) = 4
ICOUP(KM) = ICOUP(KM1)+IY(IT2+J2(KM),3)
ICOUP1(KM) = ICOUP1(KM1)+IY(IT1+J2(KM),3)
GO TO 20
17 KM = KM-1
if (KM == L) GO TO 10
GO TO 12
20 if (KM /= LN+1) GO TO 11
ITAI(ICOUP(LN+1)) = ICOUP1(LN+1)
GO TO 12
10 continue

return

end subroutine TAIL
