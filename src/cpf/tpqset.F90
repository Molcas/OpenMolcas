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
!               1986, Margareta R. A. Blomberg                         *
!***********************************************************************

subroutine TPQSET(ICASE,TPQ,IP)

implicit real*8(A-H,O-Z)
dimension TPQ(*), IOCR(100)
dimension ICASE(*)
#include "SysDef.fh"
#include "cpfmcpf.fh"
#include "spin_cpf.fh"
! Statement function
!PAM97      EXTERNAL UNPACK
!PAM97      INTEGER UNPACK
!RL   JO(L)=IAND(ISHFT(QOCC((L+29)/30),-2*((L+29)/30*30-L)),3)
!PAM97      JO(L)=UNPACK(QOCC((L+29)/30), 2*L-(2*L-1)/60*60, 2)
JO(L) = ICUNP(ICASE,L)

D0 = 0.0d0
D1 = 1.0d0
D2 = 2.0d0

IOR = 0
II1 = (IREF0-1)*LN
do I=1,LN
  JOJ = JO(II1+I)
  IOR = IOR+1
  IOCR(IOR) = JOJ
end do

IINT = IRC(4)
do IQ=1,IINT
  TPQ(IQ) = D1
  if (INCPF == 1) TPQ(IQ) = D2/N
  if ((IQ == IREF0) .or. (IP == IREF0)) TPQ(IQ) = D0
end do
if ((ISDCI == 1) .or. (INCPF == 1) .or. (IP == IREF0)) return

II = 0
IJ = 0
do I=1,LN
  JJ = (IP-1)*LN+I
  if ((JO(JJ) == IOCR(I)) .or. (JO(JJ) == 3)) GO TO 15
  if (LWSP .and. (JO(JJ)*IOCR(I) == 2)) GO TO 15
  if (II /= 0) GO TO 16
  II = I
16 IJ = I
15 continue
end do
NI = IOCR(II)
if (NI > 1) NI = NI-1
NJ = IOCR(IJ)
if (NJ > 1) NJ = NJ-1
do IQ=1,IINT
  IK = 0
  IL = 0
  do I=1,LN
    JJ = (IQ-1)*LN+I
    if ((JO(JJ) == IOCR(I)) .or. (JO(JJ) == 3)) GO TO 25
    if (LWSP .and. (JO(JJ)*IOCR(I) == 2)) GO TO 25
    if (IK /= 0) GO TO 26
    IK = I
26  IL = I
25  continue
  end do
  DIK = D0
  DIL = D0
  DJK = D0
  DJL = D0
  if (II == IK) DIK = D1
  if (II == IL) DIL = D1
  if (IJ == IK) DJK = D1
  if (IJ == IL) DJL = D1
  TPQ(IQ) = (DIK+DIL)/(D2*NI)+(DJK+DJL)/(D2*NJ)
  if (IQ == IREF0) TPQ(IQ) = D0
end do
if (IPRINT < 15) return
if (IPRINT > 5) write(6,11) (TPQ(IQ),IQ=1,IINT)
11 format(5X,'TPQ',10F5.2)

return

end subroutine TPQSET
