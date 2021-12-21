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

subroutine INT3(I,J,K,L,IT1,IT2,II,IID,JJ,JJD,JTYP,ITAI,L0,L1,L2,L3)
! K < I < J < L

use Definitions, only: wp, iwp

implicit none
integer(kind=iwp) :: I, J, K, L, IT1, IT2, II, IID, JJ, JJD, JTYP, ITAI(*), L0(*), L1(*), L2(*), L3(*)
#include "real_guga.fh"
#include "integ.fh"
integer(kind=iwp) :: ISTOP, ITAIL, ITURN, ITYP, KM, LJ, LJM, LJS

ITYP = 0
LJS = IJ(L+1)+1
LJM = IJ(L)
do LJ=LJS,LJM
  ITAIL = IX(IT2+LJ)
  if (IT1 /= IT2) call TAIL(L,LJ,ITAI,ITAIL,L0,L1,L2,L3,IT1,IT2)
  IWAY(L) = 1
32 KM = L
  J2(KM+1) = LJ
  J1(KM+1) = LJ
  call LOOP1(KM,ISTOP,IT1,IT2)
  if (ISTOP == 1) GO TO 10
41 KM = KM-1
  IWAY(KM) = 1
  if (KM == J) GO TO 51
42 call LOOP5(KM,ISTOP,IT1,IT2)
  if (ISTOP == 0) GO TO 41
  KM = KM+1
  if (KM == L) GO TO 32
  GO TO 42
51 ITURN = 0
55 IWAY(J) = 1
52 KM = J
  JM(KM) = IVF0+1
  JM1(KM) = IVF0+1
  if (ITURN == 0) call LOOP10(KM,ISTOP,IT1,IT2)
  if (ITURN == 1) call LOOP11(KM,ISTOP,IT1,IT2)
  if (ISTOP == 0) GO TO 53
  if (ITURN == 1) GO TO 54
  ITURN = 1
  GO TO 55
54 KM = KM+1
  if (KM == L) GO TO 32
  GO TO 42
53 KM = KM-1
  IWAY(KM) = 1
  if (KM == I) GO TO 61
62 JM(KM) = IVF0+1
  JM1(KM) = IVF0+1
  if (ITURN == 0) call LOOP17(KM,ISTOP,IT1,IT2)
  if (ITURN == 1) call LOOP21(KM,ISTOP,IT1,IT2)
  if (ISTOP == 0) GO TO 53
  KM = KM+1
  if (KM == J) GO TO 52
  GO TO 62
61 IWAY(I) = 1
72 KM = I
  if (ITURN == 0) call LOOP15(KM,ISTOP,IT1,IT2)
  if (ITURN == 1) call LOOP19(KM,ISTOP,IT1,IT2)
  if (ISTOP == 0) GO TO 63
  KM = KM+1
  if (KM == J) GO TO 52
  GO TO 62
63 KM = KM-1
  if (KM == 0) GO TO 20
  IWAY(KM) = 1
  if (KM == K) GO TO 74
82 call LOOP5(KM,ISTOP,IT1,IT2)
  if (ISTOP == 0) GO TO 63
  KM = KM+1
  if (KM == I) GO TO 72
  GO TO 82
74 IWAY(K) = 1
71 KM = K
  call LOOP3(KM,ISTOP,IT1,IT2)
  if (ISTOP == 1) GO TO 73
  if (abs(COUP(K)) < 1.0e-6_wp) GO TO 71
  call COMP(K,LJ,ITYP,K,IT1,IT2)
  GO TO 71
73 KM = KM+1
  if (KM == I) GO TO 72
  GO TO 82
20 call COMP1(LJ,ITYP,L,IT2,II,IID,JJ,JJD,JTYP,ITAI)
  KM = 1
  if (KM == I) GO TO 72
  GO TO 82
10 continue
end do

return

end subroutine INT3
