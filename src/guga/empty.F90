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

subroutine EMPTY(BUF,IBUF,LASTAD,SO,KBUF,NTPB)

use Definitions, only: wp, iwp

implicit none
integer(kind=iwp) :: KBUF, IBUF(KBUF+2), LASTAD(*), NTPB
real(kind=wp) :: BUF(KBUF), SO(*)
#include "SysDef.fh"
#include "real_guga.fh"
#include "integ.fh"
#include "files_addr.fh"
#include "d.fh"
integer(kind=iwp) :: I, IADR, ICLR, II, IIQQ, IJJ, IKK, IN, IND, IOFF, IQ, ISUM, ITYP, IVL, IVL0, J, JJ, JJ1, JJ2, KK, LENGTH, NBX
integer(kind=iwp), external :: ICUNP
! statement function
integer(kind=iwp) :: JO, L
JO(L) = ICUNP(ICASE,L)

ISUM = JRC(ILIM)
IOUT = 0
NMAT = 0
ICLR = NTPB
IN = ICLR+1
NBX = 0
IOFF = 0
do II=1,ISUM
  if (II > JRC(1)) GO TO 11
  IVL = IV0
  KK = II
  GO TO 15
11 if (II > JRC(2)) GO TO 12
  IVL = IV1
  KK = II-JRC(1)
  GO TO 15
12 if (II > JRC(3)) GO TO 13
  IVL = IV2
  KK = II-JRC(2)
  GO TO 15
13 IVL = IV3
  KK = II-JRC(3)
15 IOUT = IOUT+1
  ICOP1(IOUT) = 0
  if (IOUT < NBUF) GO TO 460
  ICOP1(nCOP+1) = NBUF
  call dDAFILE(Lu_10,1,COP,NCOP,IADD10)
  call iDAFILE(Lu_10,1,iCOP1,NCOP+1,IADD10)
  NMAT = NMAT+NBUF
  IOUT = 0
460 IVL0 = IV0-IVL
  !IND = KK+2**16*IVL0
  IND = ior(KK,ishft(IVL0,16))
  IOUT = IOUT+1
  ICOP1(IOUT) = IND
  if (IOUT < NBUF) GO TO 16
  ICOP1(nCOP+1) = NBUF
  call dDAFILE(Lu_10,1,COP,NCOP,IADD10)
  call iDAFILE(Lu_10,1,iCOP1,NCOP+1,IADD10)
  NMAT = NMAT+NBUF
  IOUT = 0
16 IJJ = 0
  JJ = (II-1)*LN
  do I=1,LN
    JJ1 = JO(JJ+I)
    do J=1,I
      IJJ = IJJ+1
      IN = IN+1
      if (IN <= ICLR) GO TO 100
      IN = 1
      do IIQQ=1,ICLR
        SO(IIQQ) = D0
      end do
      NBX = NBX+1
      IADR = LASTAD(NBX)
110   continue
      if (IADR == -1) GO TO 120
      ! FPS
      call dDAFILE(Lu_11,2,BUF,KBUF,IADR)
      call iDAFILE(Lu_11,2,IBUF,KBUF+2,IADR)
      LENGTH = IBUF(KBUF+1)
      IADR = IBUF(KBUF+2)
      if (LENGTH == 0) GO TO 110
      do IIQQ=1,LENGTH
        IQ = IBUF(IIQQ)-IOFF
        SO(IQ) = BUF(IIQQ)
      end do
      GO TO 110
120   IOFF = IOFF+ICLR
100   if (JJ1 == 0) GO TO 25
      JJ2 = JO(JJ+J)
      if (JJ2 == 0) GO TO 25
      ITYP = 0
      if (I == J) ITYP = 1
      IKK = IJJ
      if (ITYP == 1) IKK = I
      if (SO(IN) == D0) GO TO 25
      IOUT = IOUT+1
      !PAM96 ICOP1(IOUT) = ior(ITYP,ishft(IKK,1))
      !ICOP1(IOUT) = ITYP+2*IKK
      ICOP1(IOUT) = ior(ITYP,ishft(IKK,1))
      COP(IOUT) = SO(IN)
      if (IOUT < NBUF) GO TO 25
      ICOP1(nCOP+1) = NBUF
      call dDAFILE(Lu_10,1,COP,NCOP,IADD10)
      call iDAFILE(Lu_10,1,iCOP1,NCOP+1,IADD10)
      NMAT = NMAT+NBUF
      IOUT = 0
25    continue
    end do
  end do
end do

return

end subroutine EMPTY
