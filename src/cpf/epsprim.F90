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

subroutine EPSPRIM(JSY,INDEX,C,S,EPP)

implicit real*8(A-H,O-Z)
#include "SysDef.fh"
#include "cpfmcpf.fh"
dimension JSY(*), index(*), C(*), S(*), EPP(*)
! Statement function
!PAM97      EXTERNAL UNPACK
!PAM97      INTEGER UNPACK
!RL   JSYM(L)=IAND(ISHFT(JSY((L+19)/20),-3*((L+19)/20*20-L)),7)+1
!PAM96      JSYM(L)=UNPACK(JSY((L+9)/10),3*MOD(L-1,10)+1,3)+1
JSYM(L) = JSUNP_CPF(JSY,L)

! VALENCE

IP = IRC(1)
do I=1,IP
  EPP(I) = EPP(I)+C(I)*S(I)
end do

! SINGLES

IP = IRC(2)-IRC(1)
IN = IRC(1)
do I=1,IP
  NS1 = JSYM(IN+I)
  NSIL = MUL(NS1,LSYM)
  INUM = NVIR(NSIL)
  IST = index(IN+I)+1
  !RL call DOTPR(C(IST),1,S(IST),1,T,INUM)
  T = DDOT_(INUM,C(IST),1,S(IST),1)
  EPP(IN+I) = EPP(IN+I)+T
end do

! DOUBLES

IP = IRC(4)-IRC(2)
IN = IRC(2)
do I=1,IP
  NS1 = JSYM(IN+I)
  NSIL = MUL(NS1,LSYM)
  INUM = NNS(NSIL)
  IST = index(IN+I)+1
  !RL call DOTPR(C(IST),1,S(IST),1,T,INUM)
  T = DDOT_(INUM,C(IST),1,S(IST),1)
  EPP(IN+I) = EPP(IN+I)+T
end do

IP = IRC(4)
if (IPRINT > 5) write(6,998) (EPP(I),I=1,IP)
998 format(6X,'EPP ',5F10.6)

return

end subroutine EPSPRIM
