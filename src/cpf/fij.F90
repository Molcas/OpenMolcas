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
subroutine FIJ(ICASE,JSY,INDEX,C,S,FC,A,B,FK,DBK,ENP,EPP)

implicit real*8(A-H,O-Z)
#include "SysDef.fh"
#include "cpfmcpf.fh"
#include "files_cpf.fh"
dimension JSY(*), index(*), C(*), S(*), FC(*), A(*), B(*), FK(*), DBK(*), ENP(*), EPP(*)
dimension ICASE(*)
parameter(IPOW6=2**6,IPOW10=2**10,IPOW19=2**19)
! Statement function
JSYM(L) = JSUNP_CPF(JSY,L)

IK = 0 ! dummy initialize
NOB2 = IROW(NORBT+1)
!if (IDENS == 1) write(6,876) (FC(I),I=1,NOB2)
!876 format(1X,'FIJ',5F12.6)
ICHK = 0
if (IDENS == 1) GO TO 105
NOB2 = IROW(NORBT+1)
call SETZ(FC,NOB2)
IADD25 = 0
call dDAFILE(Lu_25,2,FC,NOB2,IADD25)
if (ITER == 1) GO TO 200
105 IADD10 = IAD10(8)
100 call dDAFILE(Lu_CIGuga,2,COP,nCOP,IADD10)
call iDAFILE(Lu_CIGuga,2,iCOP1,nCOP+1,IADD10)
LEN = ICOP1(nCOP+1)
if (LEN == 0) GO TO 100
if (LEN < 0) GO TO 200
do IN=1,LEN
  IND = ICOP1(IN)
  if (ICHK /= 0) GO TO 460
  if (IND /= 0) GO TO 11
  ICHK = 1
  GO TO 10
460 ICHK = 0
  INDI = IND
  !NI = mod(INDI,1024)
  !NK = mod(INDI/IPOW10,1024)
  NI = ibits(INDI,0,10)
  NK = ibits(INDI,10,10)
  IK = IROW(NK)+NI
  GO TO 10
11 continue
  !IVL = mod(IND,64)
  !IC2 = mod(IND/IPOW6,8192)
  !IC1 = mod(IND/IPOW19,8192)
  IVL = ibits(IND,0,6)
  IC2 = ibits(IND,6,13)
  IC1 = ibits(IND,19,13)
  COPI = COP(IN)*FC(IK)
  if (IVL /= IV0) GO TO 13
  if (IC1 /= IREF0) GO TO 16
  if (IDENS == 1) GO TO 18
  COPI = COPI/sqrt(ENP(IC2))
  S(IC2) = S(IC2)+COPI
  if (ITER == 1) GO TO 10
  EPP(IC2) = EPP(IC2)+COPI*C(IC2)
  GO TO 10
18 FC(IK) = FC(IK)+COP(IN)*C(IC1)*C(IC2)/ENP(IC2)
  GO TO 10
16 if (IC2 /= IREF0) GO TO 17
  if (IDENS == 1) GO TO 19
  COPI = COPI/sqrt(ENP(IC1))
  S(IC1) = S(IC1)+COPI
  if (ITER == 1) GO TO 10
  EPP(IC1) = EPP(IC1)+COPI*C(IC1)
  GO TO 10
19 FC(IK) = FC(IK)+COP(IN)*C(IC1)*C(IC2)/ENP(IC1)
  GO TO 10
17 if (IDENS == 1) GO TO 21
  FACS = D1
  S(IC1) = S(IC1)+FACS*COPI*C(IC2)
  S(IC2) = S(IC2)+FACS*COPI*C(IC1)
  GO TO 10
21 FC(IK) = FC(IK)+COP(IN)*C(IC1)*C(IC2)/(sqrt(ENP(IC1))*sqrt(ENP(IC2)))
  GO TO 10
13 INDA = IRC(IVL)+IC1
  INDB = IRC(IVL)+IC2
  NA = index(INDA)
  NB = index(INDB)
  NS1 = JSYM(INDA)
  NS1L = MUL(NS1,LSYM)
  INUM = NVIR(NS1L)
  if (IVL >= 2) INUM = NNS(NS1L)
  if (IDENS == 1) GO TO 15
  FACS = D1
  call DAXPY_(INUM,COPI*FACS,C(NB+1),1,S(NA+1),1)
  call DAXPY_(INUM,COPI*FACS,C(NA+1),1,S(NB+1),1)
  GO TO 10
15 TERM = DDOT_(INUM,C(NA+1),1,C(NB+1),1)
  FC(IK) = FC(IK)+COP(IN)*TERM/(sqrt(ENP(INDA))*sqrt(ENP(INDB)))
10 continue
end do
GO TO 100
!200 IF (IDENS == 1) write(6,876) (FC(I),I=1,NOB2)
200 call dAI_CPF(C)
if (ITER == 1) return
call AB(ICASE,JSY,INDEX,C,S,FC,A,B,FK,ENP)

return

! This is to allow type punning without an explicit interface
contains
subroutine dAI_CPF(C)
  use iso_c_binding
  real*8, target :: C(*)
  integer, pointer :: iC(:)
  call c_f_pointer(c_loc(C(1)),iC,[1])
  call AI_CPF(JSY,INDEX,C,S,FC,C,iC,A,B,FK,DBK,ENP,EPP,0)
  nullify(iC)
end subroutine dAI_CPF

end subroutine FIJ
