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

subroutine EPSBIS(JSY,INDX,C,W,EPB)

use Definitions, only: wp, iwp, u6, r8

implicit none
integer(kind=iwp) :: JSY(*), INDX(*)
real(kind=wp) :: C(*), W(*), EPB(*)
#include "cpfmcpf.fh"
integer(kind=iwp) :: I, IIN, INUM, IP, IST, NS1, NSIL
integer(kind=iwp), external :: JSUNP_CPF
real(kind=r8), external :: DDOT_
! Statement function
!PAM97      EXTERNAL UNPACK
!PAM97      INTEGER UNPACK
!RL   JSYM(L)=IAND(ISHFT(JSY((L+19)/20),-3*((L+19)/20*20-L)),7)+1
!PAM96      JSYM(L)=UNPACK(JSY((L+9)/10),3*MOD(L-1,10)+1,3)+1
integer(kind=iwp) :: JSYM, L
JSYM(L) = JSUNP_CPF(JSY,L)

call SETZ(EPB,IRC(4))
if ((ICPF == 1) .or. (ISDCI == 1) .or. (INCPF == 1)) return

! VALENCE

IP = IRC(1)
do I=1,IP
  EPB(I) = C(I)*W(I)
end do

! SINGLES

IP = IRC(2)-IRC(1)
IIN = IRC(1)
do I=1,IP
  !FUE IND = IND+1
  NS1 = JSYM(IIN+I)
  NSIL = MUL(NS1,LSYM)
  INUM = NVIR(NSIL)
  IST = INDX(IIN+I)+1
  !RL call DOTPR(C(IST),1,W(IST),1,EPB(IIN+I),INUM)
  EPB(IIN+I) = DDOT_(INUM,C(IST),1,W(IST),1)
end do

! DOUBLES

IP = IRC(4)-IRC(2)
IIN = IRC(2)
do I=1,IP
  NS1 = JSYM(IIN+I)
  NSIL = MUL(NS1,LSYM)
  INUM = NNS(NSIL)
  IST = INDX(IIN+I)+1
  !RL call DOTPR(C(IST),1,W(IST),1,EPB(IIN+I),INUM)
  EPB(IIN+I) = DDOT_(INUM,C(IST),1,W(IST),1)
end do

IP = IRC(4)
if (IPRINT > 5) write(u6,998) (EPB(I),I=1,IP)

return

998 format(6X,'EPB ',5F10.6)

end subroutine EPSBIS
