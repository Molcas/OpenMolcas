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

subroutine TAB2(NREF,IOCR,nIOCR,L0,L1,L2,L3,INTNUM,LV,LSYM,ICIALL,IFCORE,ICOR,NONE,JONE)

use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp) :: NREF, nIOCR, IOCR(nIOCR), L0(*), L1(*), L2(*), L3(*), INTNUM, LV, LSYM, ICIALL, IFCORE, ICOR(*), NONE, JONE(*)
#include "real_guga.fh"
#include "integ.fh"
integer(kind=iwp) :: I, I0, I1, IA1, IAC, IAT, IB1, IBMAX, IBS, IBT, IEL, II, IIJ, IIJF, IIM, IIM2, IJD, IJFL, IJFS, IJL, IJR, &
                     IJRL, IJS, IJS1, IL, IN, INUM, IORB(MXVERT), ISTA, ISTOP, ISUM, ITTT, IUT, IUT1, J, J11, J3, J4, JJ, JJ1, &
                     JJ2, JL, JMAX, K, K00(MXVERT), K11(MXVERT), K22(MXVERT), K33(MXVERT), L00(MXVERT), L11(MXVERT), L22(MXVERT), &
                     L33(MXVERT), LN1, NAC, NACU, NIJ, NIJ1
real(kind=wp) :: FB, FBB

IEL = 2
if (IFIRST /= 0) IEL = 1

! NUMBER OF ACTIVE ELECTRONS
NAC = N-2*NIORB

! UPPER LIMIT FOR NUMBER OF ELECTRONS IN ACTIVE SPACE
NACU = NAC+IEL

IUT = 0
IB(1) = int(2*S)
IA(1) = int(N-2*S)/2
IJ(LN+1) = 0
IJ(LN) = 1
NIJ = 1
IJR = 1
IJS = 2
IJRL = IJR
IORB(1) = 0
do II=1,LN
  IIM = LN-II+1-LV
  IAC = IIM-NIORB-1
  IIM2 = (IIM-1)*2

  ! S=0

16 INUM = N-2*IA(IJR)-IB(IJR)
  if ((INTNUM == 0) .or. (IIM <= 0)) GO TO 109
  if (IAC >= 0) GO TO 409
  if (IB(IJR) == 0) GO TO 109
  if (INUM+IIM2 == N-2) GO TO 11
  GO TO 109
409 if (IB(IJR) > IAC+3) GO TO 11
109 NIJ = NIJ+1
  IA(NIJ) = IA(IJR)
  IB(NIJ) = IB(IJR)
  IORB(NIJ) = IORB(IJR)+2
  if (IIM <= 0) GO TO 11
  call CHEL(IA(NIJ),IB(NIJ),IIM,IEL,ISTOP)
  if (ISTOP == 1) NIJ = NIJ-1

  ! S=1

11 if (IB(IJR) == 0) GO TO 12
  if ((INTNUM == 0) .or. (IIM <= 0)) GO TO 112
  if (IAC < 0) GO TO 112
  if (INUM+1 > NACU) GO TO 12
  if (INUM+1 /= NACU) GO TO 112
  IBS = IB(IJR)-1
  if ((IBS /= 0) .and. (IBS /= 2)) GO TO 12
112 NIJ = NIJ+1
  IA(NIJ) = IA(IJR)
  IB(NIJ) = IB(IJR)-1
  IORB(NIJ) = IORB(IJR)+1
  if (IIM <= 0) GO TO 12
  call CHEL(IA(NIJ),IB(NIJ),IIM,IEL,ISTOP)
  if (ISTOP == 1) NIJ = NIJ-1

  ! S=2

12 if (IA(IJR) == 0) GO TO 13
  if ((INTNUM == 0) .or. (IIM <= 0)) GO TO 113
  if (IAC < 0) GO TO 413
  if (IB(IJR)+1 > IAC+3) GO TO 13
  if (INUM+1 > NACU) GO TO 13
  if (INUM+1 /= NACU) GO TO 113
  IBS = IB(IJR)+1
  if ((IBS /= 0) .and. (IBS /= 2)) GO TO 13
  GO TO 113
413 if (IB(IJR) >= 2) GO TO 13
113 NIJ = NIJ+1
  IA(NIJ) = IA(IJR)-1
  IB(NIJ) = IB(IJR)+1
  IORB(NIJ) = IORB(IJR)+1
  if (IIM <= 0) GO TO 13
  call CHEL(IA(NIJ),IB(NIJ),IIM,IEL,ISTOP)
  if (ISTOP == 1) NIJ = NIJ-1

  ! S=3

13 if (IA(IJR) == 0) GO TO 14
  if (INTNUM == 0) GO TO 114
  if (IAC < 0) GO TO 114
  if (IB(IJR) > IAC+3) GO TO 14
  if (INUM+2 > NACU) GO TO 14
  if (INUM+2 /= NACU) GO TO 114
  IBS = IB(IJR)
  if ((IBS /= 0) .and. (IBS /= 2)) GO TO 14
114 NIJ = NIJ+1
  IA(NIJ) = IA(IJR)-1
  IB(NIJ) = IB(IJR)
  IORB(NIJ) = IORB(IJR)
  if (IIM <= 0) GO TO 14
  call CHEL(IA(NIJ),IB(NIJ),IIM,IEL,ISTOP)
  if (ISTOP == 1) NIJ = NIJ-1

14 if (IJR == IJRL) GO TO 15
  IJR = IJR+1
  GO TO 16
  ! DELETE VERTICES
15 NIJ1 = NIJ-1
  IN = IJS
  IUT = IJS
  if (NIJ1 < IJS) GO TO 21
  do IJD=IJS,NIJ1
    JJ1 = NIJ-IJD+IJS-1
    J = JJ1+1
    do K=IJS,JJ1
      if (IA(J) /= IA(K)) GO TO 25
      if (IB(J) /= IB(K)) GO TO 25
      GO TO 26
25    continue
    end do
    GO TO 20
26  IA(J) = -1
    IB(J) = -1
20  continue
  end do
  ! PACK VERTICES
  IJS1 = IJS+1
  do J=IJS1,NIJ
    if (IA(J) /= -1) GO TO 31
    if (IB(J) /= -1) GO TO 31
    IN = IN+1
    GO TO 30
31  IN = IN+1
    IUT = IUT+1
    IA(IUT) = IA(IN)
    IB(IUT) = IB(IN)
    IORB(IUT) = IORB(IN)
30  continue
  end do
  ! ORDER VERTICES
  IUT1 = IUT-1
  if (IUT1 < IJS) GO TO 21
  do J=IJS,IUT1
    J11 = J+1
    do K=J11,IUT
      if (IA(J)-IA(K) < 0) then
        GO TO 43
      else if (IA(J)-IA(K) == 0) then
        GO TO 44
      else
        GO TO 42
      end if
44    if (IB(J) > IB(K)) GO TO 42
43    IAT = IA(J)
      IBT = IB(J)
      IA(J) = IA(K)
      IB(J) = IB(K)
      IA(K) = IAT
      IB(K) = IBT
42    continue
    end do
  end do
21 if (II /= LN) IJ(LN-II) = IUT
  IJR = IJS
  IJS = IUT+1
  IJRL = IUT
  NIJ = IUT
end do
if (N /= 2) GO TO 35
IUT = IUT+1
IA(IUT) = 0
IB(IUT) = 0
IA(IUT-1) = 0
IB(IUT-1) = 1
IA(IUT-2) = 0
IB(IUT-2) = 2
35 JJ2 = 0
do II=1,LN
  I = LN-II+1-LV
  IIM2 = (I-1)*2
  I0 = I+LV
  JJ1 = IJ(I0+1)+1
  JJ2 = IJ(I0)
  J3 = JJ2+1
  if (I0 /= 1) J4 = IJ(I0-1)
  if (I0 == 1) J4 = IUT
  ! DETERMINE CASE DOWN
  do J=JJ1,JJ2
    IA1 = IA(J)
    IB1 = IB(J)
    INUM = 2*IA1+IB1
    K00(J) = 0
    K11(J) = 0
    K22(J) = 0
    K33(J) = 0
    do JJ=J3,J4
      if (IA1 == IA(JJ)) GO TO 61
      if ((IA1-IA(JJ)) /= 1) GO TO 60
      if (IB1 == IB(JJ)) GO TO 62
      if ((IB(JJ)-IB1) /= 1) GO TO 60
      if ((I > NIORB) .or. (INTNUM == 0)) GO TO 59
      if (I <= 0) GO TO 59
      if (IB1 >= 2) GO TO 60
59    K22(J) = JJ
      GO TO 60
62    K33(J) = JJ
      GO TO 60
61    if (IB1 == IB(JJ)) GO TO 63
      if ((IB1-IB(JJ)) /= 1) GO TO 60
      K11(J) = JJ
      GO TO 60
63    if ((I > NIORB) .or. (INTNUM == 0)) GO TO 64
      if (I <= 0) GO TO 64
      if (IB1 == 0) GO TO 64
      if (INUM-IIM2 == 2) GO TO 60
      if (IA1 < I-1) GO TO 60
64    K00(J) = JJ
60    continue
    end do
  end do
  ! DETERMINE CASE UP
  do J=J3,J4
    IA1 = IA(J)
    IB1 = IB(J)
    INUM = 2*IA1+IB1
    L00(J) = 0
    L11(J) = 0
    L22(J) = 0
    L33(J) = 0
    do JJ=JJ1,JJ2
      if (IA(JJ) == IA1) GO TO 81
      if ((IA(JJ)-IA1) /= 1) GO TO 80
      if (IB(JJ) == IB1) GO TO 82
      if ((IB1-IB(JJ)) /= 1) GO TO 80
      if ((I > NIORB) .or. (INTNUM == 0)) GO TO 79
      if (I <= 0) GO TO 79
      if (IB1 >= 3) GO TO 80
79    L22(J) = JJ
      GO TO 80
82    L33(J) = JJ
      GO TO 80
81    if (IB(JJ) == IB1) GO TO 83
      if ((IB(JJ)-IB1) /= 1) GO TO 80
      L11(J) = JJ
      GO TO 80
83    if ((I > NIORB) .or. (INTNUM == 0)) GO TO 84
      if (I <= 0) GO TO 84
      if (IB1 == 0) GO TO 84
      if (INUM-IIM2 == 2) GO TO 80
      if (IA1 < I-1) GO TO 80
84    L00(J) = JJ
80    continue
    end do
  end do
end do
IV0 = IUT
IV1 = IUT-1
IV2 = IUT-2
IV3 = IUT-3
K00(IUT) = 0
K11(IUT) = 0
K22(IUT) = 0
K33(IUT) = 0
K00(IUT+1) = 0
K11(IUT+1) = 0
K22(IUT+1) = 0
K33(IUT+1) = 0
if (ICIALL /= 0) call CIALL(LSYM,NREF,IOCR,nIOCR,L00,L11,L22,L33,LV)
call DELTAB(NREF,IOCR,L0,L1,L2,L3,INTNUM,LV,IFCORE,ICOR,NONE,JONE,K00,K11,K22,K33,L00,L11,L22,L33)
do I=1,ILIM
  ISTA = (I-1)*MXVERT
  if (IPRINT >= 5) write(IW,101)
101 format(///,6X,'TAB2',//,8X,'J',8X,'A',3X,'B',7X,'K0',2X,'K1',2X,'K2',2X,'K3',2X,'L0',2X,'L1',2X,'L2',2X,'L3',/)
  if (IPRINT >= 5) write(IW,100) (J,IA(J),IB(J),K0(ISTA+J),K1(ISTA+J),K2(ISTA+J),K3(ISTA+J),L0(ISTA+J),L1(ISTA+J),L2(ISTA+J), &
                                  L3(ISTA+J),J=1,IUT)
100 format(6X,I3,5X,2I4,5X,8I4)
end do
IBMAX = 0
do J=1,IUT
  if (IB(J) > IBMAX) IBMAX = IB(J)
end do
IUT1 = IUT-1
do IL=1,ILIM
  ISTA = (IL-1)*MXVERT
  IX(ISTA+1) = 1
  do II=1,LN
    if (II == LN) GO TO 420
    I = LN-II
    IJL = IJ(I)
    IJS = IJ(I+1)+1
    GO TO 430
420 IJL = IUT
    IJS = IUT-3
    if (IFIRST /= 0) IJS = IUT-1
430 do J=IJS,IJL
      ISUM = 0
      if (L0(ISTA+J) == 0) GO TO 441
      ISUM = ISUM+IX(ISTA+L0(ISTA+J))
441   if (L1(ISTA+J) == 0) GO TO 442
      IY(ISTA+L1(ISTA+J),1) = ISUM
      ISUM = ISUM+IX(ISTA+L1(ISTA+J))
442   if (L2(ISTA+J) == 0) GO TO 443
      IY(ISTA+L2(ISTA+J),2) = ISUM
      ISUM = ISUM+IX(ISTA+L2(ISTA+J))
443   if (L3(ISTA+J) == 0) GO TO 444
      IY(ISTA+L3(ISTA+J),3) = ISUM
      ISUM = ISUM+IX(ISTA+L3(ISTA+J))
444   IX(ISTA+J) = ISUM
    end do
  end do
end do
do I=1,ILIM
  ISTA = (I-1)*MXVERT
  if (IPRINT >= 5) write(IW,102)
102 format(///,6X,'INDEX TABLE',//,8X,'J',8X,'Y1',3X,'Y2',3X,'Y3',9X,'X',/)
  if (IPRINT >= 5) write(IW,200) (J,IY(ISTA+J,1),IY(ISTA+J,2),IY(ISTA+J,3),IX(ISTA+J),J=1,IUT)
200 format(6X,I3,5X,3I5,5X,I5)
end do
if (IPRINT >= 2) write(IW,210) IUT
210 format(/,6X,'NUMBER OF VERTICES',I10)
write(IW,214)
214 format(///,6X,'INTERNAL CONFIGURATIONS (FORMAL)')

if (IFIRST == 0) then
! This is a normal calculation, both singles and doubles included.
  write(IW,215) (IX(IUT+1-ITTT+MXVERT*(ITTT-1)),ITTT=1,4)
215 format(/,6X,'NUMBER OF VALENCE STATES',I16,/,6X,'NUMBER OF DOUBLET COUPLED SINGLES',I7,/,6X, &
           'NUMBER OF TRIPLET COUPLED DOUBLES',I7,/,6X,'NUMBER OF SINGLET COUPLED DOUBLES',I7)
else
  ! "FIRST" keyword has been given. Then this is just a so-called
  ! first-order CI, i.e., singles only.
  write(IW,215) (IX(IUT+1-ITTT+MXVERT*(ITTT-1)),ITTT=1,2)
  !216   FORMAT(/,6X,'NUMBER OF VALENCE STATES',I16,/,6X,'NUMBER OF DOUBLET COUPLED SINGLES',I7)
end if

IRC(1) = IX(IUT)
do I=2,ILIM
  ISTA = (I-1)*MXVERT
  IRC(I) = IX(ISTA+IUT+1-I)+IRC(I-1)
end do
ISUM = IRC(ILIM)
if (ISUM > LIX) then
  write(u6,*) 'Tab2: ISUM > LIX'
  write(u6,*) 'ISUM,LIX=',ISUM,LIX
  call Abend()
end if
do I=1,10
  FBB = I-1
  FB = FBB/(FBB+1)
  BS1(I) = sqrt(FB)
  FB = (FBB+2)/(FBB+1)
  BS2(I) = sqrt(FB)
  if (I > 1) BS3(I) = D1/BS1(I)
  BS4(I) = D1/BS2(I)
  FB = FBB*FBB-1
  if (I > 1) BL1(I) = sqrt(FB)/FBB
  FB = (FBB+2)**2-1
  BL2(I) = sqrt(FB)/(FBB+2)
end do
! PUT ZEROS IN VECTORS
do I=1,LN
  COUP(I) = D0
  COUP1(I) = D0
end do
IN = 0
do I=1,LN
  II = LN-I+1
  IJS = IJ(II+1)+1
  IJL = IJ(II)
  IJFS = IJF(II+1)+1
  IJFL = IJF(II)
  do IIJ=IJS,IJL
    do IIJF=IJFS,IJFL
      if (IA(IIJ) /= IAF(IIJF)) GO TO 330
      if (IB(IIJ) /= IBF(IIJF)) GO TO 330
      IPO(IIJ) = IIJF
      GO TO 320
330   continue
    end do
320 continue
  end do
end do
IPO(IUT) = IJF(1)+1
if (IPRINT >= 10) write(IW,350) (IPO(J),J=1,IUT)
350 format(/,6X,5I5)
JMAX = 0
LN1 = LN+1
do I=2,LN1
  I1 = I-1
  JL = IJ(I1)-IJ(I)
  if (JL < JMAX) GO TO 400
  JMAX = JL
400 continue
end do

if (IPRINT >= 2) then
  write(IW,411) JMAX
  write(IW,412) IBMAX,MAXB
411 format(/,6X,'NUMBER OF VERTICES IN ONE ROW',I6,6X,'(PRESENT LIMIT 31)')
  if (JMAX > 31) call Abend()
412 format(6X,'MAXIMUM B VALUE',I20,6X,'(PRESENT LIMIT ',I5,')')
end if
if (IBMAX > MAXB) call Abend()

return

end subroutine TAB2
