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

subroutine AIBJ(L0,L1,L2,L3,ITAI)

implicit real*8(A-H,O-Z)
#include "SysDef.fh"
dimension L0(*), L1(*), L2(*), L3(*), ITAI(*)
#include "real_guga.fh"
#include "integ.fh"
#include "files_addr.fh"
#include "d.fh"
dimension NUMM(7)
! statement function
JO(L) = ICUNP(ICASE,L)

IC1 = 0    ! dummy initialize
IC2 = 0    ! dummy initialize
COPLA = D0 ! dummy initialize
do I=1,7
  NUMM(I) = 0
end do
IOUT = 0
NMAT = 0
do NI=1,LN
  do NJ=1,NI
    I = ICH(NI)
    J = ICH(NJ)
    if (I > J) GO TO 19
    I = ICH(NJ)
    J = ICH(NI)
19  LTYP = 0
    IOUT = IOUT+1
    ICOP1(IOUT) = 0
    if (IOUT < NBUF) GO TO 460
    ICOP1(nCOP+1) = NBUF
    call dDAFILE(Lu_10,1,COP,NCOP,IADD10)
    call iDAFILE(Lu_10,1,iCOP1,NCOP+1,IADD10)
    NMAT = NMAT+NBUF
    IOUT = 0
460 IOUT = IOUT+1
    !PAM96 ICOP1(IOUT) = ior(I,ishft(J,10))
    !ICOP1(IOUT) = I+2**10*J
    ICOP1(IOUT) = ior(I,ishft(J,10))
    if (IOUT < NBUF) GO TO 11
    ICOP1(nCOP+1) = NBUF
    call dDAFILE(Lu_10,1,COP,NCOP,IADD10)
    call iDAFILE(Lu_10,1,iCOP1,NCOP+1,IADD10)
    NMAT = NMAT+NBUF
    IOUT = 0
11  IJS = IJ(I+1)+1
    IJM = IJ(I)
    ! FIRST ORDER INTERACTION
    ! TRIPLET-VALENCE INTERACTIONS
    JTURN = 1
    ITURN = 0
    ITT1 = 2
    ITT2 = 0
150 IT1 = ITT1*MXVERT
    !ulf IT2 = ITT2*300
    IT2 = ITT2*MXVERT
    II = 0
    IID = 0
    if (ITT2 == 0) GO TO 149
    II = IRC(ITT2)
    IID = JRC(ITT2)
149 JJ = IRC(ITT1)
    JJD = JRC(ITT1)
    ITYP = JTURN
    if (ITURN /= 0) ITYP = JTURN-2
    do IJJ=IJS,IJM
      ITAIL = IX(IT2+IJJ)
      if (IT1 /= IT2) call TAIL(I,IJJ,ITAI,ITAIL,L0,L1,L2,L3,IT1,IT2)
      IWAY(I) = 1
32    KM = I
      J2(KM+1) = IJJ
      J1(KM+1) = IJJ
      if (I == J) GO TO 51
      call LOOP1(KM,ISTOP,IT1,IT2)
      if (ISTOP == 1) GO TO 30
41    KM = KM-1
      IWAY(KM) = 1
      if (KM == J) GO TO 51
42    call LOOP5(KM,ISTOP,IT1,IT2)
      if (ISTOP == 0) GO TO 41
      KM = KM+1
      if (KM == I) GO TO 32
      GO TO 42
51    IWAY(J) = 1
52    KM = J
      JM(KM) = IVF0+1
      JM1(KM) = IVF0+1
      IABIJ = 0
      if (I == J) GO TO 12
      if (ITURN == 0) call LOOP10(KM,ISTOP,IT1,IT2)
      IFAI = 1
      if (ITURN == 1) call LOOP13(KM,ISTOP,IFAI,IT1,IT2)
      if (ISTOP == 1) GO TO 14
      if (J1(KM) /= J2(KM)) GO TO 13
      IABIJ = 1
      IC1 = ICOUP(KM)
      IC2 = ICOUP1(KM)
      KM1 = KM+1
      IDIF = IA(J1(KM1))-IA(J2(KM1))
      if (IWAY(J) == 2) COPLA = COUP(J+1)
      if ((IDIF == 0) .and. (IWAY(J) == 5)) COPLA = COUP(J+1)*BS4(IB(J2(J+1))+1)
      if ((IDIF == 1) .and. (IWAY(J) == 4)) COPLA = COUP(J+1)*BS3(IB(J2(J+1))+1)
      GO TO 53
12    if (ITURN == 0) call LOOP7(KM,ISTOP,IT1,IT2)
      IFAI = 1
      if (ITURN == 1) call LOOP26(KM,ISTOP,IFAI,IT1,IT2)
      LTYP = 1
      if (IWAY(I) == 5) LTYP = 0
13    if (ISTOP == 0) GO TO 53
14    if (I == J) GO TO 30
      KM = KM+1
      if (KM == I) GO TO 32
      GO TO 42
53    KM = KM-1
      if (KM == 0) GO TO 61
      IWAY(KM) = 1
62    JM(KM) = IVF0+1
      JM1(KM) = IVF0+1
      if (ITURN == 0) call LOOP17(KM,ISTOP,IT1,IT2)
      IFAI = 0
      if (J1(KM+1) /= J2(KM+1)) GO TO 165
      if (I == J) GO TO 162
      if (IABIJ == 0) GO TO 165
      IC11 = ICOUP(KM+1)-IC1
      IC22 = ICOUP1(KM+1)-IC2
      if (IC11 == IC22) IFAI = 1
      GO TO 165
162   if (ICOUP(KM+1) /= ICOUP1(KM+1)) GO TO 165
      IFAI = 1
165   if (ITURN == 1) call LOOP23(KM,ISTOP,IFAI,IT1,IT2)
      if (ISTOP == 0) GO TO 53
      KM = KM+1
      if (KM == J) GO TO 52
      GO TO 62
61    COPL = COUP(1)
      if ((JTURN >= 6) .or. (JTURN == 3)) COPL = COUP1(1)
      IFAB = 0
      if (JTURN /= 5) GO TO 63
      if ((I == J) .and. (LTYP == 1)) GO TO 72
63    if (ITT1 /= ITT2) GO TO 70
      if (I /= J) GO TO 65
      if ((LTYP == 1) .and. (ICOUP(1) > ICOUP1(1))) GO TO 72
65    if (IABIJ == 0) GO TO 70
      IC11 = ICOUP(1)-IC1
      IC22 = ICOUP1(1)-IC2
      if (IC11 == IC22) IFAB = 1
70    COPLA0 = COPLA
      if (IFAB == 0) COPLA0 = D0
      do IN=1,ITAIL
        ICP1 = ICOUP(1)+IN
        JND1 = JNDX(II+ICP1)
        if (JND1 == 0) GO TO 80
        ICP1 = JND1-IID
        if (ITT1 /= ITT2) GO TO 289
        IN2 = IN
        GO TO 288
289     IN2 = ITAI(IN)
        if (IN2 == 0) GO TO 80
288     ICP2 = ICOUP1(1)+IN2
        JND2 = JNDX(JJ+ICP2)
        if (JND2 == 0) GO TO 80
        ICP2 = JND2-JJD
        if ((ITT1 /= ITT2) .or. (ICP1 /= ICP2)) GO TO 100
        JJ1 = (JJD+ICP1-1)*LN+I
        JOJ = JO(JJ1)
        if (JOJ > 1) JOJ = JOJ-1
        COPLA0 = JOJ-2
        IFAB = 1
100     IOUT = IOUT+1
        NUMM(JTURN) = NUMM(JTURN)+1
        COP(IOUT) = COPL
        !PAM96 IND1 = ior(IFAB,ishft(ITURN,1))
        !PAM96 IND2 = ior(IND1,ishft(ITYP,2))
        !PAM96 IND3 = ior(IND2,ishft(ICP1,5))
        !PAM96 ICOP1(IOUT) = ior(IND3,ishft(ICP2,18))
        !IND1 = IFAB+2**1*ITURN
        !IND2 = IND1+2**2*ITYP
        !IND3 = IND2+2**5*ICP1
        !ICOP1(IOUT) = IND3+2**18*ICP2
        IND1 = ior(IFAB,ishft(ITURN,1))
        IND2 = ior(IND1,ishft(ITYP,2))
        IND3 = ior(IND2,ishft(ICP1,5))
        ICOP1(IOUT) = ior(IND3,ishft(ICP2,18))
        if (IOUT < NBUF) GO TO 71
        ICOP1(nCOP+1) = NBUF
        call dDAFILE(Lu_10,1,COP,NCOP,IADD10)
        call iDAFILE(Lu_10,1,iCOP1,NCOP+1,IADD10)
        NMAT = NMAT+NBUF
        IOUT = 0
71      if (IFAB == 0) GO TO 80
        IOUT = IOUT+1
        NUMM(JTURN) = NUMM(JTURN)+1
        COP(IOUT) = COPLA0
        ICOP1(IOUT) = 1
        if (IOUT < NBUF) GO TO 80
        ICOP1(nCOP+1) = NBUF
        call dDAFILE(Lu_10,1,COP,NCOP,IADD10)
        call iDAFILE(Lu_10,1,iCOP1,NCOP+1,IADD10)
        NMAT = NMAT+NBUF
        IOUT = 0
80      continue
      end do
72    KM = KM+1
      if (KM == J) GO TO 52
      GO TO 62
30    continue
    end do
    GO TO(151,152,153,154,155,156,20),JTURN
    ! SINGLET-VALENCE INTERACTIONS
151 JTURN = 2
    ITT1 = 3
    GO TO 150
    ! TRIPLET-TRIPLET INTERACTIONS
152 JTURN = 3
    ITURN = 1
    ITT1 = 2
    ITT2 = 2
    GO TO 150
    ! SINGLET-SINGLET INTERACTIONS
153 JTURN = 4
    ITT1 = 3
    ITT2 = 3
    GO TO 150
    ! TRIPLET-SINGLET INTERACTIONS
154 JTURN = 5
    ITT1 = 2
    ITT2 = 3
    GO TO 150
    ! SINGLET-TRIPLET INTERACTIONS
155 JTURN = 6
    ITT1 = 3
    ITT2 = 2
    GO TO 150
    ! DOUBLET-DOUBLET INTERACTIONS
156 JTURN = 7
    ITT1 = 1
    ITT2 = 1
    GO TO 150
20  continue
  end do
end do
ICOP1(nCOP+1) = IOUT
call dDAFILE(Lu_10,1,COP,NCOP,IADD10)
call iDAFILE(Lu_10,1,iCOP1,NCOP+1,IADD10)
NMAT = NMAT+IOUT
ICOP1(nCOP+1) = -1
call dDAFILE(Lu_10,1,COP,NCOP,IADD10)
call iDAFILE(Lu_10,1,iCOP1,NCOP+1,IADD10)
write(IW,600) NMAT
600 format(/,6X,'COEFFICIENTS FOR AIBJ',I9)
write(IW,610) (NUMM(I),I=1,7)
610 format(6X,'DIFFERENT TYPES',7I9)

return

end subroutine AIBJ
