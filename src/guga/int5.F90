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

subroutine INT5(I,J,L)
! I < J < L   I == K
implicit real*8(A-H,O-Z)
#include "real_guga.fh"
#include "integ.fh"

ITYP = 0
LJS = IJ(L+1)+1
LJM = IJ(L)
do ITT=1,ILIM
  IT1 = (ITT-1)*MXVERT
  IT2 = IT1
  do LJ=LJS,LJM
    IWAY(L) = 1
32  KM = L
    J2(KM+1) = LJ
    J1(KM+1) = LJ
    call LOOP1(KM,ISTOP,IT1,IT2)
    if (ISTOP == 1) GO TO 10
41  KM = KM-1
    IWAY(KM) = 1
    if (KM == J) GO TO 51
42  call LOOP5(KM,ISTOP,IT1,IT2)
    if (ISTOP == 0) GO TO 41
    KM = KM+1
    if (KM == L) GO TO 32
    GO TO 42
51  ITURN = 0
55  IWAY(J) = 1
52  KM = J
    JM(KM) = IVF0+1
    JM1(KM) = IVF0+1
    if (ITURN == 0) call LOOP10(KM,ISTOP,IT1,IT2)
    IFAI = 0
    if (ITURN == 1) call LOOP13(KM,ISTOP,IFAI,IT1,IT2)
    if (ISTOP == 0) GO TO 53
    if (ITURN == 1) GO TO 54
    ITURN = 1
    GO TO 55
54  KM = KM+1
    if (KM == L) GO TO 32
    GO TO 42
53  KM = KM-1
    IWAY(KM) = 1
    if (KM == I) GO TO 71
62  JM(KM) = IVF0+1
    JM1(KM) = IVF0+1
    if (ITURN == 0) call LOOP17(KM,ISTOP,IT1,IT2)
    IFAI = 0
    if (ITURN == 1) call LOOP23(KM,ISTOP,IFAI,IT1,IT2)
    if (ISTOP == 0) GO TO 53
    KM = KM+1
    if (KM == J) GO TO 52
    GO TO 62
71  KM = I
    if (ITURN == 0) call LOOP14(KM,ISTOP,IT1,IT2)
    if (ITURN == 1) call LOOP22(KM,ISTOP,IT1,IT2)
    if (ISTOP == 1) GO TO 73
    if (abs(COUP(I)) < 1.D-06) GO TO 71
    call COMP(I,LJ,ITYP,I,IT1,IT2)
    GO TO 71
73  KM = KM+1
    if (KM == J) GO TO 52
    GO TO 62
10  continue
  end do
end do

return

end subroutine INT5
