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

subroutine INT8(I,J,L,IT1,IT2,II,IID,JJT,JJD,JTYP,ITAI,L0,L1,L2,L3)
! K == L , I < J

use Definitions, only: wp, iwp

#include "intent.fh"

implicit none
integer(kind=iwp), intent(in) :: I, J, L, IT1, IT2, II, IID, JJT, JJD, JTYP, L0(*), L1(*), L2(*), L3(*)
integer(kind=iwp), intent(_OUT_) :: ITAI(*)
#include "real_guga.fh"
#include "integ.fh"
integer(kind=iwp) :: ISTOP, ITAIL, ITYP, JJ, JJM, JJS, KM
real(kind=wp) :: FAC

ITYP = 0
if ((L < I) .or. (L > J)) ITYP = 1
if (I == 0) ITYP = 0
if ((I == 0) .and. (L > J)) ITYP = 1
FAC = D1
JJS = IJ(J+1)+1
JJM = IJ(J)
do JJ=JJS,JJM
  ITAIL = IX(IT2+JJ)
  if (IT1 /= IT2) call TAIL(J,JJ,ITAI,ITAIL,L0,L1,L2,L3,IT1,IT2)
  IWAY(J) = 1
32 KM = J
  J2(KM+1) = JJ
  J1(KM+1) = JJ
  call LOOP1(KM,ISTOP,IT1,IT2)
  if (ISTOP == 1) GO TO 10
41 KM = KM-1
  if (KM == 0) GO TO 20
  IWAY(KM) = 1
  if (KM == I) GO TO 51
42 call LOOP5(KM,ISTOP,IT1,IT2)
  if (ISTOP == 1) GO TO 43
  if (KM /= L) GO TO 41
  if (IWAY(KM) == 2) GO TO 42
  FAC = D1
  if (IWAY(KM) == 5) FAC = D2
  GO TO 41
43 KM = KM+1
  if (KM == J) GO TO 32
  GO TO 42
51 IWAY(I) = 1
52 KM = I
  call LOOP3(KM,ISTOP,IT1,IT2)
  if (ISTOP == 1) GO TO 53
  if (abs(COUP(I)) < 1.0e-6_wp) GO TO 52
  COUP(I) = FAC*COUP(I)
  call COMP(I,JJ,ITYP,L,IT1,IT2)
  GO TO 52
53 KM = KM+1
  if (KM == J) GO TO 32
  GO TO 42
20 COUP(1) = FAC*COUP(1)
  call COMP1(JJ,ITYP,L,IT2,II,IID,JJT,JJD,JTYP,ITAI)
  KM = 1
  if (KM == J) GO TO 32
  GO TO 42
10 continue
end do

return

end subroutine INT8
