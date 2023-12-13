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

use guga_global, only: IADD10, ICASE, ILIM, IOUT, IV0, JRC, LN, Lu_10, Lu_11, NBUF, NMAT
use guga_util_global, only: COP, ICOP1, nCOP
use Constants, only: Zero
use Definitions, only: wp, iwp

#include "intent.fh"

implicit none
integer(kind=iwp), intent(in) :: KBUF, LASTAD(*), NTPB
real(kind=wp), intent(out) :: BUF(KBUF)
integer(kind=iwp), intent(out) :: IBUF(KBUF+2)
real(kind=wp), intent(_OUT_) :: SO(*)
integer(kind=iwp) :: I, IADR, ICLR, II, IIQQ, IJJ, IKK, IN_, IND, IOFF, IQ, ISUM, ITYP, IVL, IVL0, J, JJ, JJ1, JJ2, KK, LENGTH, NBX
integer(kind=iwp), external :: ICUNP

ISUM = JRC(ILIM)
IOUT = 0
NMAT = 0
ICLR = NTPB
IN_ = ICLR+1
NBX = 0
IOFF = 0
do II=1,ISUM
  if (II <= JRC(1)) then
    IVL = IV0
    KK = II
  else if (II <= JRC(2)) then
    IVL = IV0-1
    KK = II-JRC(1)
  else if (II <= JRC(3)) then
    IVL = IV0-2
    KK = II-JRC(2)
  else
    IVL = IV0-3
    KK = II-JRC(3)
  end if
  IOUT = IOUT+1
  ICOP1(IOUT) = 0
  if (IOUT >= NBUF) then
    ICOP1(nCOP+1) = NBUF
    call dDAFILE(Lu_10,1,COP,NCOP,IADD10)
    call iDAFILE(Lu_10,1,iCOP1,NCOP+1,IADD10)
    NMAT = NMAT+NBUF
    IOUT = 0
  end if
  IVL0 = IV0-IVL
  !IND = KK+2**16*IVL0
  IND = ior(KK,ishft(IVL0,16))
  IOUT = IOUT+1
  ICOP1(IOUT) = IND
  if (IOUT >= NBUF) then
    ICOP1(nCOP+1) = NBUF
    call dDAFILE(Lu_10,1,COP,NCOP,IADD10)
    call iDAFILE(Lu_10,1,iCOP1,NCOP+1,IADD10)
    NMAT = NMAT+NBUF
    IOUT = 0
  end if
  IJJ = 0
  JJ = (II-1)*LN
  do I=1,LN
    JJ1 = ICUNP(ICASE,JJ+I)
    do J=1,I
      IJJ = IJJ+1
      IN_ = IN_+1
      if (IN_ > ICLR) then
        IN_ = 1
        SO(1:ICLR) = Zero
        NBX = NBX+1
        IADR = LASTAD(NBX)
        do while (IADR /= -1)
          ! FPS
          call dDAFILE(Lu_11,2,BUF,KBUF,IADR)
          call iDAFILE(Lu_11,2,IBUF,KBUF+2,IADR)
          LENGTH = IBUF(KBUF+1)
          IADR = IBUF(KBUF+2)
          do IIQQ=1,LENGTH
            IQ = IBUF(IIQQ)-IOFF
            SO(IQ) = BUF(IIQQ)
          end do
        end do
        IOFF = IOFF+ICLR
      end if
      if (JJ1 == 0) cycle
      JJ2 = ICUNP(ICASE,JJ+J)
      if (JJ2 == 0) cycle
      ITYP = 0
      if (I == J) ITYP = 1
      IKK = IJJ
      if (ITYP == 1) IKK = I
      if (SO(IN_) == Zero) cycle
      IOUT = IOUT+1
      !PAM96 ICOP1(IOUT) = ior(ITYP,ishft(IKK,1))
      !ICOP1(IOUT) = ITYP+2*IKK
      ICOP1(IOUT) = ior(ITYP,ishft(IKK,1))
      COP(IOUT) = SO(IN_)
      if (IOUT < NBUF) cycle
      ICOP1(nCOP+1) = NBUF
      call dDAFILE(Lu_10,1,COP,NCOP,IADD10)
      call iDAFILE(Lu_10,1,iCOP1,NCOP+1,IADD10)
      NMAT = NMAT+NBUF
      IOUT = 0
    end do
  end do
end do

return

end subroutine EMPTY
