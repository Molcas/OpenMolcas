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

subroutine DELTAB(NREF,IOCR,L0,L1,L2,L3,INTNUM,LV,IFCORE,ICOR,NONE,JONE,K00,K11,K22,K33,L00,L11,L22,L33)

use Definitions, only: iwp, u6

implicit none
integer(kind=iwp) :: NREF, IOCR(*), L0(*), L1(*), L2(*), L3(*), INTNUM, LV, IFCORE, ICOR(*), NONE, JONE(*), K00(*), K11(*), &
                     K22(*), K33(*), L00(*), L11(*), L22(*), L33(*)
#include "real_guga.fh"
#include "integ.fh"
integer(kind=iwp) :: I, IBS, IDIF, IEL, IIJ, IJJ, INHOLE, IOC(55), IPART, IREF, IRR, ISP(55), ISTA, JHOLE, JJ1, JPART, K, &
                     K0M(MXVERT), K1M(MXVERT), K2M(MXVERT), K3M(MXVERT), KM, KM1, L0M(MXVERT), L1M(MXVERT), L2M(MXVERT), &
                     L3M(MXVERT), LNS, LSYM, NCORR, NSJ

do I=1,4*MXVERT
  K0(I) = 0
  K1(I) = 0
  K2(I) = 0
  K3(I) = 0
  L0(I) = 0
  L1(I) = 0
  L2(I) = 0
  L3(I) = 0
end do
IBS = 0
IEL = 2
LSYM = 1
LNS = NIORB+LV+1
IRR = 0
do I=LNS,LN
  IRR = IRR+1
  if (IOCR(IRR) == 1) LSYM = MUL(LSYM,NSM(I))
end do
do IIJ=1,ILIM
  ISTA = (IIJ-1)*MXVERT
  do I=1,IV0
    K0M(I) = 0
    K1M(I) = 0
    K2M(I) = 0
    K3M(I) = 0
    L0M(I) = 0
    L1M(I) = 0
    L2M(I) = 0
    L3M(I) = 0
  end do
  IJJ = IV0+1-IIJ
  KM = 1
  J2(KM) = IJJ
11 KM = KM+1
  IWAY(KM) = 0
12 KM1 = KM-1
  if ((L00(J2(KM1)) == 0) .or. (IWAY(KM) >= 1)) GO TO 14
  J2(KM) = L00(J2(KM1))
  IWAY(KM) = 1
  IOC(KM1) = 0
  ISP(KM1) = 0
  GO TO 20
14 if ((L11(J2(KM1)) == 0) .or. (IWAY(KM) >= 2)) GO TO 15
  J2(KM) = L11(J2(KM1))
  IWAY(KM) = 2
  IOC(KM1) = 1
  ISP(KM1) = 1
  GO TO 20
15 if ((L22(J2(KM1)) == 0) .or. (IWAY(KM) >= 3)) GO TO 16
  J2(KM) = L22(J2(KM1))
  IWAY(KM) = 3
  IOC(KM1) = 1
  ISP(KM1) = 2
  GO TO 20
16 if ((L33(J2(KM1)) == 0) .or. (IWAY(KM) >= 4)) GO TO 17
  J2(KM) = L33(J2(KM1))
  IWAY(KM) = 4
  IOC(KM1) = 2
  ISP(KM1) = 3
  GO TO 20
17 KM = KM-1
  if (KM == 1) GO TO 210
  GO TO 12
20 if (KM1 == NIORB+LV) IBS = IB(J2(KM))
  if (KM /= LN+1) GO TO 11
  NSJ = 1
  INHOLE = 0
  do I=1,LN
    if (IOC(I) == 1) NSJ = MUL(NSJ,NSM(I))
    if ((I <= NIORB+LV) .and. (I > LV)) INHOLE = INHOLE+2-IOC(I)
  end do
  ! STRIKE OUT INTERNAL CONFIGURATIONS
  IPART = 0
  if (IIJ > 1) IPART = IPART+1
  if (IIJ > 2) IPART = IPART+1
  JJ1 = 0
  do IREF=1,NREF
    JHOLE = 0
    JPART = IPART
    do I=1,LN
      if (I > LV) GO TO 250
      IDIF = IOC(I)
      GO TO 251
250   if (I > NIORB+LV) GO TO 252
      IDIF = IOC(I)-2
      GO TO 251
252   JJ1 = JJ1+1
      if (IOC(I) == IOCR(JJ1)) GO TO 112
      IDIF = IOC(I)-IOCR(JJ1)
251   if (IDIF > 0) GO TO 114
      JHOLE = JHOLE-IDIF
      GO TO 112
114   JPART = JPART+IDIF
112   continue
    end do
    if (JPART /= JHOLE) then
      write(u6,*) 'DeltaB: JPART.NE.JHOLE'
      write(u6,*) 'JPART,JHOLE=',JPART,JHOLE
      write(u6,*) 'iREF=',iREF
      call Abend()
    end if
    if (JPART <= IEL) GO TO 113
  end do
  GO TO 12
113 if ((IPART == 0) .and. (NSJ /= LSYM)) GO TO 12
  if ((IPART /= 2) .or. (INTNUM == 0)) GO TO 115
  ! INTERACTING SPACE
  if ((INHOLE == 2) .and. (IBS /= 0)) GO TO 12
  ! NO CORE-CORE CORRELATION
115 if (IFCORE == 0) GO TO 116
  NCORR = 0
  do I=1,LN
    if (ICOR(I) == 0) GO TO 117
    NCORR = NCORR+2-IOC(I)
117 continue
  end do
  if (NCORR > 1) GO TO 12
  ! SINGLY OCCUPIED ORBITALS
116 if (NONE == 0) GO TO 118
  do I=1,NONE
    if (IOC(JONE(I)) /= 1) GO TO 12
  end do
118 do K=1,LN
    if (ISP(K)-1 < 0) then
      GO TO 131
    else if (ISP(K)-1 == 0) then
      GO TO 132
    else
      GO TO 133
    end if
131 K0M(J2(K+1)) = K00(J2(K+1))
    L0M(J2(K)) = L00(J2(K))
    GO TO 130
132 K1M(J2(K+1)) = K11(J2(K+1))
    L1M(J2(K)) = L11(J2(K))
    GO TO 130
133 if (ISP(K) == 3) GO TO 134
    K2M(J2(K+1)) = K22(J2(K+1))
    L2M(J2(K)) = L22(J2(K))
    GO TO 130
134 K3M(J2(K+1)) = K33(J2(K+1))
    L3M(J2(K)) = L33(J2(K))
130 continue
  end do
  GO TO 12
210 do I=1,IV0
    K0(ISTA+I) = K0M(I)
    K1(ISTA+I) = K1M(I)
    K2(ISTA+I) = K2M(I)
    K3(ISTA+I) = K3M(I)
    L0(ISTA+I) = L0M(I)
    L1(ISTA+I) = L1M(I)
    L2(ISTA+I) = L2M(I)
    L3(ISTA+I) = L3M(I)
  end do
end do

return

end subroutine DELTAB
