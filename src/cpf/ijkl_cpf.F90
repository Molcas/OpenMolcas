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

!pgi$g opt=1
subroutine IJKL_CPF(JSY,INDEX,C,S,FIJKL,BUFIN,IBUFIN,ENP,EPP)

implicit real*8(A-H,O-Z)
#include "SysDef.fh"
#include "cpfmcpf.fh"
#include "files_cpf.fh"
dimension JSY(*), index(*), C(*), S(*), FIJKL(*), BUFIN(*), IBUFIN(*), ENP(*), EPP(*)
parameter(IPOW8=2**8,IPOW16=2**16,IPOW24=2**24)
parameter(IPOW6=2**6,IPOW13=2**13,IPOW19=2**19)
! Statement function
JSYM(L) = JSUNP_CPF(JSY,L)

FINI = 0.0d0 ! dummy initialize
NCONF = JSC(4)
ICHK = 0
NIJ = IROW(LN+1)
NIJKL = NIJ*(NIJ+1)/2
do I=1,NIJKL
  FIJKL(I) = D0
end do
KKBUF0 = (RTOI*(KBUFF1+2)-2)/(RTOI+1)
KKBUF1 = RTOI*KKBUF0+KKBUF0+1
KKBUF2 = KKBUF1+1
IADR = LASTAD(1)
201 call iDAFILE(Lu_TiABCI,2,IBUFIN,KKBUF2,IADR)
LENGTH = IBUFIN(KKBUF1)
IADR = IBUFIN(KKBUF2)
if (LENGTH == 0) GO TO 209
call SCATTER(LENGTH,FIJKL,IBUFIN(RTOI*KKBUF0+1),BUFIN)
209 if (IADR /= -1) GO TO 201
IADD10 = IAD10(5)
100 call dDAFILE(Lu_CIGuga,2,COP,nCOP,IADD10)
call iDAFILE(Lu_CIGuga,2,iCOP1,nCOP+1,IADD10)
LEN = ICOP1(nCOP+1)
if (LEN == 0) GO TO 100
if (LEN < 0) GO TO 200
do IN=1,LEN
  IND = ICOP1(IN)
  if (ICHK /= 0) GO TO 460
  if (IND /= 0) GO TO 22
  ICHK = 1
  GO TO 10
460 ICHK = 0
  INDI = IND
  !IP = mod(INDI,IPOW8)
  !JP = mod(INDI/IPOW8,IPOW8)
  !KP = mod(INDI/IPOW16,IPOW8)
  !LP = mod(INDI/IPOW24,IPOW8)
  IP = ibits(INDI,0,8)
  JP = ibits(INDI,8,8)
  KP = ibits(INDI,16,8)
  LP = ibits(INDI,24,8)
  NIJ = IROW(IP)+JP
  NKL = IROW(KP)+LP
  IND = NIJ*(NIJ-1)/2+NKL
  FINI = FIJKL(IND)
  GO TO 10
22 if (abs(FINI) < 1.d-06) GO TO 10
  !PAM97 IVL = iand(IND,63)
  !PAM97 IC2 = iand(ishft(IND,-6),8191)
  !PAM97 IC1 = iand(ishft(IND,-19),8191)
  !IVL = mod(IND,IPOW6)
  !IC2 = mod(IND/IPOW6,IPOW13)
  !IC1 = mod(IND/IPOW19,IPOW13)
  IVL = ibits(IND,0,6)
  IC2 = ibits(IND,6,13)
  IC1 = ibits(IND,19,13)
  COPI = COP(IN)*FINI
  if (IVL /= 0) GO TO 13
  if (IC1 /= IREF0) GO TO 16
  COPI = COPI/sqrt(ENP(IC2))
  S(IC2) = S(IC2)+COPI
  if (ITER == 1) GO TO 10
  EPP(IC2) = EPP(IC2)+COPI*C(IC2)
  GO TO 10
16 if (IC2 /= IREF0) GO TO 17
  COPI = COPI/sqrt(ENP(IC1))
  S(IC1) = S(IC1)+COPI
  if (ITER == 1) GO TO 10
  EPP(IC1) = EPP(IC1)+COPI*C(IC1)
  GO TO 10
17 FACS = D1
  S(IC1) = S(IC1)+FACS*COPI*C(IC2)
  S(IC2) = S(IC2)+FACS*COPI*C(IC1)
  GO TO 10
13 INDA = IRC(IVL)+IC1
  INDB = IRC(IVL)+IC2
  FACS = D1
  NA = index(INDA)
  NB = index(INDB)
  NS1 = JSYM(INDA)
  NS1L = MUL(NS1,LSYM)
  INUM = NVIR(NS1L)
  if (IVL >= 2) INUM = NNS(NS1L)
  call DAXPY_(INUM,COPI*FACS,C(NB+1),1,S(NA+1),1)
  call DAXPY_(INUM,COPI*FACS,C(NA+1),1,S(NB+1),1)
10 continue
end do
GO TO 100

200 return

end subroutine IJKL_CPF
