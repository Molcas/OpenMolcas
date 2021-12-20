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

subroutine INT1(I,J,K,L,IT1,IT2,II,IID,JJ,JJD,JTYP,ITAI,L0,L1,L2,L3)

implicit real*8(A-H,O-Z)
dimension ITAI(*), L0(*), L1(*), L2(*), L3(*)
! I < J < K < L
#include "real_guga.fh"
#include "integ.fh"

LJS = IJ(L+1)+1
LJM = IJ(L)
ITYP = 0
do LJ=LJS,LJM
  ITAIL = IX(IT2+LJ)
  if (IT1 /= IT2) call TAIL(L,LJ,ITAI,ITAIL,L0,L1,L2,L3,IT1,IT2)
  IWAY(L) = 1
32 KM = L
  J2(KM+1) = LJ
  J1(KM+1) = LJ
  call LOOP1(L,ISTOP,IT1,IT2)
  if (ISTOP == 1) GO TO 10
41 KM = KM-1
  IWAY(KM) = 1
  if (KM == K) GO TO 51
42 call LOOP5(KM,ISTOP,IT1,IT2)
  if (ISTOP == 0) GO TO 41
  KM = KM+1
  if (KM == L) GO TO 32
  GO TO 42
51 IWAY(K) = 1
52 KM = K
  call LOOP3(KM,ISTOP,IT1,IT2)
  if (ISTOP == 0) GO TO 30
  KM = KM+1
  if (KM == L) GO TO 32
  GO TO 42
30 KM = KM-1
  IWAY(KM) = 1
  if (KM == J) GO TO 133
62 call PATH(KM,ISTOP,IT1,IT2)
  if (ISTOP == 0) GO TO 30
  KM = KM+1
  if (KM == K) GO TO 52
  GO TO 62
133 IWAY(J) = 1
132 KM = J
  if (IT1 >= IT2) call LOOP1(KM,ISTOP,IT1,IT2)
  if (IT2 > IT1) call LOOP2(KM,ISTOP,IT1,IT2)
  if (ISTOP == 0) GO TO 141
  KM = KM+1
  if (KM == K) GO TO 52
  GO TO 62
141 KM = KM-1
  if (KM == 0) GO TO 20
  IWAY(KM) = 1
  if (KM == I) GO TO 151
142 if (IT1 >= IT2) call LOOP5(KM,ISTOP,IT1,IT2)
  if (IT2 > IT1) call LOOP6(KM,ISTOP,IT1,IT2)
  if (ISTOP == 0) GO TO 141
  KM = KM+1
  if (KM == J) GO TO 132
  GO TO 142
151 IWAY(I) = 1
152 KM = I
  if (IT1 >= IT2) call LOOP3(KM,ISTOP,IT1,IT2)
  if (IT2 > IT1) call LOOP4(KM,ISTOP,IT1,IT2)
  if (ISTOP == 1) GO TO 153
  COUP(I) = COUP(I)*COUP(K)
  if (abs(COUP(I)) < 1.D-06) GO TO 152
  ICPI = ICOUP(I)
  ICP1I = ICOUP1(I)
  ICOUP(I) = ICOUP(J+1)+ICPI
  ICOUP1(I) = ICOUP1(J+1)+ICP1I
  call COMP(I,LJ,ITYP,I,IT1,IT2)
  ICOUP(I) = ICOUP(J+1)+ICP1I
  ICOUP1(I) = ICOUP1(J+1)+ICPI
  call COMP(I,LJ,ITYP,I,IT1,IT2)
  GO TO 152
153 KM = KM+1
  if (KM == J) GO TO 132
  GO TO 142
20 COUP(1) = COUP(1)*COUP(K)
  if (abs(COUP(1)) < 1.D-06) GO TO 154
  ICOUP(1) = ICOUP(J+1)+ICOUP(1)
  ICOUP1(1) = ICOUP1(J+1)+ICOUP1(1)
  call COMP1(LJ,ITYP,L,IT2,II,IID,JJ,JJD,JTYP,ITAI)
154 KM = 1
  if (KM == J) GO TO 132
  GO TO 142
10 continue
end do

return

end subroutine INT1
