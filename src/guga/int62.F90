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

subroutine INT62(I,K,L,IT1,IT2,II,IID,JJ,JJD,JTYP,ITAI,L0,L1,L2,L3)
! I < K < L  J == L

use Definitions, only: wp, iwp

implicit none
integer(kind=iwp) :: I, K, L, IT1, IT2, II, IID, JJ, JJD, JTYP, ITAI(*), L0(*), L1(*), L2(*), L3(*)
#include "real_guga.fh"
#include "integ.fh"
integer(kind=iwp) :: ISTOP, ITAIL, ITURN, ITYP, KM, LJ, LJM, LJS

LJS = IJ(L+1)+1
LJM = IJ(L)
do LJ=LJS,LJM
  ITAIL = IX(IT2+LJ)
  if (IT1 /= IT2) call TAIL(L,LJ,ITAI,ITAIL,L0,L1,L2,L3,IT1,IT2)
  IWAY(L) = 1
32 KM = L
  J2(KM+1) = LJ
  J1(KM+1) = LJ
  JM(KM) = IVF0+1
  JM1(KM) = IVF0+1
  call LOOP8(KM,ISTOP,IT1,IT2)
  if (ISTOP == 1) GO TO 10
  ITYP = 0
  if (IWAY(L) /= 5) ITYP = 2
53 KM = KM-1
  IWAY(KM) = 1
  if (KM == K) GO TO 61
62 JM(KM) = IVF0+1
  JM1(KM) = IVF0+1
  call LOOP21(KM,ISTOP,IT1,IT2)
  if (ISTOP == 0) GO TO 53
  KM = KM+1
  if (KM == L) GO TO 32
  GO TO 62
61 ITURN = 0
65 IWAY(K) = 1
72 KM = K
  if (ITURN == 0) call LOOP20(KM,ISTOP,IT1,IT2)
  if (ITURN == 1) call LOOP19(KM,ISTOP,IT1,IT2)
  if (ISTOP == 0) GO TO 63
  if (ITURN == 1) GO TO 64
  ITURN = 1
  GO TO 65
64 KM = KM+1
  if (KM == L) GO TO 32
  GO TO 62
63 KM = KM-1
  if (KM == 0) GO TO 20
  IWAY(KM) = 1
  if (KM == I) GO TO 74
82 if (ITURN == 0) call LOOP6(KM,ISTOP,IT1,IT2)
  if (ITURN == 1) call LOOP5(KM,ISTOP,IT1,IT2)
  if (ISTOP == 0) GO TO 63
  KM = KM+1
  if (KM == K) GO TO 72
  GO TO 82
74 IWAY(I) = 1
71 KM = I
  if (ITURN == 0) call LOOP4(KM,ISTOP,IT1,IT2)
  if (ITURN == 1) call LOOP3(KM,ISTOP,IT1,IT2)
  if (ISTOP == 1) GO TO 73
  if (abs(COUP(I)) < 1.0e-6_wp) GO TO 71
  if ((ITYP == 2) .and. (ICOUP1(I) < ICOUP(I))) GO TO 71
  call COMP(I,LJ,ITYP,I,IT1,IT2)
  GO TO 71
73 KM = KM+1
  if (KM == K) GO TO 72
  GO TO 82
20 if ((JTYP > 3) .and. (ITYP == 2)) GO TO 21
  call COMP1(LJ,ITYP,L,IT2,II,IID,JJ,JJD,JTYP,ITAI)
21 KM = 1
  if (KM == K) GO TO 72
  GO TO 82
10 continue
end do

return

end subroutine INT62
