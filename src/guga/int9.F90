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

subroutine INT9(I,J,L,IT1,IT2,II,IID,JJ,JJD,JTYP,ITAI,L0,L1,L2,L3)
! I < L . CASE 1 J=K=L , CASE 2 I=J=K .

use Definitions, only: wp, iwp

implicit none
integer(kind=iwp) :: I, J, L, IT1, IT2, II, IID, JJ, JJD, JTYP, ITAI(*), L0(*), L1(*), L2(*), L3(*)
#include "real_guga.fh"
#include "integ.fh"
integer(kind=iwp) :: ISTOP, ITAIL, ITYP, KM, LJ, LJM, LJS
real(kind=wp) :: FAC

ITYP = 0
FAC = D1
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
  if (I == J) GO TO 41
  FAC = D1
  if ((IWAY(KM) == 2) .or. (IWAY(KM) == 4)) FAC = D0
41 KM = KM-1
  if (KM == 0) GO TO 20
  IWAY(KM) = 1
  if (KM == I) GO TO 51
42 call LOOP5(KM,ISTOP,IT1,IT2)
  if (ISTOP == 0) GO TO 41
  KM = KM+1
  if (KM == L) GO TO 32
  GO TO 42
51 IWAY(I) = 1
52 KM = I
  call LOOP3(KM,ISTOP,IT1,IT2)
  if (ISTOP == 1) GO TO 53
  if (I /= J) GO TO 54
  FAC = D1
  if (IWAY(KM) == 2) FAC = D0
54 if (FAC == D0) GO TO 52
  if (abs(COUP(I)) < 1.0e-6_wp) GO TO 52
  call COMP(I,LJ,ITYP,I,IT1,IT2)
  GO TO 52
53 KM = KM+1
  if (KM == L) GO TO 32
  GO TO 42
20 if (FAC == D0) GO TO 33
  call COMP1(LJ,ITYP,L,IT2,II,IID,JJ,JJD,JTYP,ITAI)
33 KM = 1
  if (KM == L) GO TO 32
  GO TO 42
10 continue
end do

return

end subroutine INT9
