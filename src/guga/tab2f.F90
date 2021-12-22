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

subroutine TAB2F(IVER,LV)

use Definitions, only: iwp, u6

implicit none
integer(kind=iwp), intent(in) :: IVER, LV
#include "real_guga.fh"
#include "integ.fh"
integer(kind=iwp) :: I, IA1, IAT, IB1, IBT, IEL, II, IIM, IJD, IJR, IJRL, IJS, IJS1, IN_, IORB(MXVERT), ISTOP, IUT, IUT1, J, J11, &
                     J3, J4, JJ, JJ1, JJ2, K, NIJ, NIJ1, nijj

nijj = 0
IEL = 2
if (IFIRST /= 0) IEL = 1
IEL = IEL+1
IUT = 0
IBF(1) = int(2*S)
IAF(1) = int(N-2*S)/2
IJF(LN+1) = 0
IJF(LN) = 1
NIJ = 1
IJR = 1
IJS = 2
IJRL = IJR
IORB(1) = 0
do II=1,LN
  IIM = LN-II+1-LV
  ! S=0
16 NIJ = NIJ+1
  !if (NIJ > IVER)
  nijj = max(nij,nijj)
  IAF(NIJ) = IAF(IJR)
  IBF(NIJ) = IBF(IJR)
  IORB(NIJ) = IORB(IJR)+2
  if (IIM <= 0) GO TO 11
  call CHEL(IAF(NIJ),IBF(NIJ),IIM,IEL,ISTOP)
  if (ISTOP == 1) NIJ = NIJ-1
11 if (IBF(IJR) == 0) GO TO 12
  ! S=1
  NIJ = NIJ+1
  !if (NIJ > IVER)
  nijj = max(nij,nijj)
  IAF(NIJ) = IAF(IJR)
  IBF(NIJ) = IBF(IJR)-1
  IORB(NIJ) = IORB(IJR)+1
  if (IIM <= 0) GO TO 12
  call CHEL(IAF(NIJ),IBF(NIJ),IIM,IEL,ISTOP)
  if (ISTOP == 1) NIJ = NIJ-1
12 if (IAF(IJR) == 0) GO TO 13
  ! S=2
  NIJ = NIJ+1
  !if (NIJ > IVER)
  nijj = max(nij,nijj)
  IAF(NIJ) = IAF(IJR)-1
  IBF(NIJ) = IBF(IJR)+1
  IORB(NIJ) = IORB(IJR)+1
  if (IIM <= 0) GO TO 13
  call CHEL(IAF(NIJ),IBF(NIJ),IIM,IEL,ISTOP)
  if (ISTOP == 1) NIJ = NIJ-1
13 if (IAF(IJR) == 0) GO TO 14
  ! S=3
  NIJ = NIJ+1
  !if (NIJ > IVER)
  nijj = max(nij,nijj)
  IAF(NIJ) = IAF(IJR)-1
  IBF(NIJ) = IBF(IJR)
  IORB(NIJ) = IORB(IJR)
  if (IIM <= 0) GO TO 14
  call CHEL(IAF(NIJ),IBF(NIJ),IIM,IEL,ISTOP)
  if (ISTOP == 1) NIJ = NIJ-1
14 if (IJR == IJRL) GO TO 15
  IJR = IJR+1
  GO TO 16
  ! DELETE VERTICES
15 continue
  NIJ1 = NIJ-1
  IN_ = IJS
  IUT = IJS
  if (NIJ1 < IJS) GO TO 21
  do IJD=IJS,NIJ1
    JJ1 = NIJ-IJD+IJS-1
    J = JJ1+1
    do K=IJS,JJ1
      if (IAF(J) /= IAF(K)) GO TO 25
      if (IBF(J) /= IBF(K)) GO TO 25
      GO TO 26
25    continue
    end do
    GO TO 20
26  IAF(J) = -1
    IBF(J) = -1
20  continue
  end do
  ! PACK VERTICES
  IJS1 = IJS+1
  do J=IJS1,NIJ
    if (IAF(J) /= -1) GO TO 31
    if (IBF(J) /= -1) GO TO 31
    IN_ = IN_+1
    GO TO 30
31  IN_ = IN_+1
    IUT = IUT+1
    IAF(IUT) = IAF(IN_)
    IBF(IUT) = IBF(IN_)
    IORB(IUT) = IORB(IN_)
30  continue
  end do
  ! ORDER VERTICES
  IUT1 = IUT-1
  if (IUT1 < IJS) GO TO 21
  do J=IJS,IUT1
    J11 = J+1
    do K=J11,IUT
      if (IAF(J)-IAF(K) < 0) then
        GO TO 43
      else if (IAF(J)-IAF(K) == 0) then
        GO TO 44
      else
        GO TO 42
      end if
44    if (IBF(J) > IBF(K)) GO TO 42
43    IAT = IAF(J)
      IBT = IBF(J)
      IAF(J) = IAF(K)
      IBF(J) = IBF(K)
      IAF(K) = IAT
      IBF(K) = IBT
42    continue
    end do
  end do
21 if (II /= LN) IJF(LN-II) = IUT
  IJR = IJS
  IJS = IUT+1
  IJRL = IUT
  NIJ = IUT
end do
JJ2 = 0
do II=1,LN
  I = LN-II+1
  JJ1 = IJF(I+1)+1
  JJ2 = IJF(I)
  J3 = JJ2+1
  if (I /= 1) J4 = IJF(I-1)
  if (I == 1) J4 = IUT
  ! DETERMINE CASE DOWN
  do J=JJ1,JJ2
    IA1 = IAF(J)
    IB1 = IBF(J)
    K0F(J) = 0
    K1F(J) = 0
    K2F(J) = 0
    K3F(J) = 0
    do JJ=J3,J4
      if (IA1 == IAF(JJ)) GO TO 61
      if ((IA1-IAF(JJ)) /= 1) GO TO 60
      if (IB1 == IBF(JJ)) GO TO 62
      if ((IBF(JJ)-IB1) /= 1) GO TO 60
      K2F(J) = JJ
      GO TO 60
62    K3F(J) = JJ
      GO TO 60
61    if (IB1 == IBF(JJ)) GO TO 63
      if ((IB1-IBF(JJ)) /= 1) GO TO 60
      K1F(J) = JJ
      GO TO 60
63    K0F(J) = JJ
60    continue
    end do
  end do
end do
IVF0 = IUT
IVF1 = IUT-1
IVF2 = IUT-2
IVF3 = IUT-3
K0F(IUT) = 0
K1F(IUT) = 0
K2F(IUT) = 0
K3F(IUT) = 0
K0F(IUT+1) = 0
K1F(IUT+1) = 0
K2F(IUT+1) = 0
K3F(IUT+1) = 0
if (IPRINT >= 5) write(IW,101)
101 format(///,6X,'TAB2F',//,8X,'J',7X,'AF',2X,'BF',6X,'K0F',1X,'K1F',1X,'K2F',1X,'K3F',/)
if (IPRINT >= 5) write(IW,100) (J,IAF(J),IBF(J),K0F(J),K1F(J),K2F(J),K3F(J),J=1,IUT)
100 format(6X,I3,5X,2I4,5X,4I4)
write(u6,*) ' Number of vertices',nijj,iut
if (nijj > iver) go to 300
return
300 write(IW,310) IVER
310 format(/,6X,'NUMBER OF VERTICES EXCEEDS',I7)
call Abend()

end subroutine TAB2F
