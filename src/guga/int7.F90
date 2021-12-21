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

subroutine INT7(I,K,L,IDIAG,BUFOUT,INDOUT,ICAD,IBUFL,KBUF,NTPB)
! I < L  I == K  J == L

use Definitions, only: wp, iwp

implicit none
integer(kind=iwp) :: I, K, L, IDIAG, INDOUT(*), ICAD(*), IBUFL(*), KBUF, NTPB
real(kind=wp) :: BUFOUT(*)
#include "SysDef.fh"
#include "real_guga.fh"
#include "integ.fh"
#include "files_addr.fh"
#include "d.fh"
integer(kind=iwp) :: IAD110, ICP, ICPP, ICQ, IDIV, IJJ, IN, IPOS, ISTOP, ISU, ISUM, IT1, IT2, ITAIL, ITT, ITURN, ITYP, IVL, JND1, &
                     KBUF0, KBUF1, KBUF2, KM, LJ, LJM, LJS, NBN

IJJ = 0 ! dummy initialize
KBUF0 = RTOI*KBUF
KBUF1 = KBUF0+KBUF+1
KBUF2 = KBUF1+1
IDIV = RTOI
ITYP = 0
if (IDIAG == 1) IJJ = L*(L-1)/2+K
LJS = IJ(L+1)+1
LJM = IJ(L)
do ITT=1,ILIM
  IT1 = (ITT-1)*MXVERT
  IT2 = IT1
  do LJ=LJS,LJM
    ITURN = 0
    if (IDIAG == 1) ITURN = 1
33  IWAY(L) = 1
32  KM = L
    J2(KM+1) = LJ
    J1(KM+1) = LJ
    JM(KM) = IVF0+1
    JM1(KM) = IVF0+1
    if (ITURN == 0) call LOOP7(KM,ISTOP,IT1,IT2)
    if (ITURN == 1) call LOOP8(KM,ISTOP,IT1,IT2)
    if (ISTOP == 1) GO TO 54
    if ((IDIAG == 1) .and. (J1(KM) /= J2(KM))) GO TO 32
    GO TO 53
54  if (ITURN == 1) GO TO 10
    ITURN = 1
    GO TO 33
53  KM = KM-1
    IWAY(KM) = 1
    if (KM == I) GO TO 71
62  JM(KM) = IVF0+1
    JM1(KM) = IVF0+1
    if (ITURN == 0) call LOOP17(KM,ISTOP,IT1,IT2)
    if (ITURN == 1) call LOOP21(KM,ISTOP,IT1,IT2)
    if (ISTOP == 1) GO TO 55
    if ((IDIAG == 1) .and. (J1(KM) /= J2(KM))) GO TO 62
    GO TO 53
55  KM = KM+1
    if (KM == L) GO TO 32
    GO TO 62
71  KM = I
    if (ITURN == 0) call LOOP14(KM,ISTOP,IT1,IT2)
    if (ITURN == 1) call LOOP18(KM,ISTOP,IT1,IT2)
    if (ISTOP == 1) GO TO 73
    if (abs(COUP(I)) < 1.0e-6_wp) GO TO 71
    if (IDIAG == 0) GO TO 105
    if (ICOUP1(I) /= ICOUP(I)) GO TO 71
    GO TO 106
105 if ((ITURN == 0) .or. (IWAY(L) == 5)) GO TO 106
    if (ICOUP1(I) < ICOUP(I)) GO TO 71
    if (ICOUP1(I) == ICOUP(I)) GO TO 71
106 if (ITURN == 0) COUP(I) = COUP(I)/D2
    if (IDIAG == 1) GO TO 25
    call COMP(I,LJ,ITYP,I,IT1,IT2)
    GO TO 71
25  KM = KM-1
    if (KM == 0) GO TO 26
    IWAY(KM) = 1
27  call PATH(KM,ISTOP,IT1,IT2)
    if (ISTOP == 0) GO TO 25
    KM = KM+1
    if (KM == I) GO TO 71
    GO TO 27
26  IVL = J2(1)
    ITAIL = IX(IT1+LJ)
    ISUM = IV0-IVL
    ISU = 0
    if (ISUM /= 0) ISU = IRC(ISUM)
    do IN=1,ITAIL
      JND1 = JNDX(ISU+ICOUP(1)+IN)
      if (JND1 == 0) GO TO 104
      IPOS = (JND1-1)*LNP+IJJ
      NBN = (IPOS-1)/NTPB+1
      IBUFL(NBN) = IBUFL(NBN)+1
      ICQ = ICAD(NBN)
      ICP = ICQ/IDIV+IBUFL(NBN)
      BUFOUT(ICP) = COUP(I)
      ICPP = ICQ+KBUF0+IBUFL(NBN)
      INDOUT(ICPP) = IPOS
      if (IBUFL(NBN) < KBUF) GO TO 104
      INDOUT(ICQ+KBUF1) = KBUF
      IAD110 = IADD11
      call iDAFILE(Lu_11,1,INDOUT(ICQ+1),KBUF2,IADD11)
      INDOUT(ICQ+KBUF2) = IAD110
      IBUFL(NBN) = 0
104   continue
    end do
    if (I == 1) GO TO 71
    KM = 1
    GO TO 27
73  KM = KM+1
    if (KM == L) GO TO 32
    GO TO 62
10  continue
  end do
end do

return

end subroutine INT7
