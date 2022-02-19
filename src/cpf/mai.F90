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

subroutine MAI(JSY,INDEX,C,S,FC,BUFIN,IBUFIN,A,B,FK,DBK,W,THET,ENP,EPP,NII,KTYP)
! KTYP=0  ,  (A/I)   INTEGRALS
! KTYP=1  ,  (AI/JK) INTEGRALS

implicit real*8(A-H,O-Z)
#include "SysDef.fh"
#include "cpfmcpf.fh"
#include "files_cpf.fh"
dimension JSY(*), index(*), C(*), S(*), FC(*), BUFIN(*), IBUFIN(*), A(*), B(*), FK(*), DBK(*), W(*), THET(NII,NII), ENP(*), EPP(*)
dimension IPOB(9)
parameter(IPOW6=2**6,IPOW13=2**13,IPOW19=2**19)
parameter(IPOW10=2**10,IPOW20=2**20)
! Statement function
JSYM(L) = JSUNP_CPF(JSY,L)

!if (IDENS == 1) write(6,876) (FC(I),I=1,NOB2)
!876 format(1X,'AI',5F12.6)
NK = 0 ! dummy initialize
NSK = 0 ! dummy initialize
INUM = IRC(4)-IRC(3)
call MPSQ2(C,S,W,MUL,INDEX,JSY,NDIAG,INUM,IRC(3),LSYM,NVIRT,SQ2)
NVT = IROW(NVIRT+1)
ICHK = 0
IJOLD = 0
NOB2 = IROW(NORBT+1)
NOT2 = IROW(LN+1)
NOTT = 2*NOT2
NOVST = LN*NVIRT+1+NVT
LBUF0 = RTOI*LBUF
LBUF1 = LBUF0+LBUF+1
LBUF2 = LBUF1+1
if (KTYP == 0) IADD10 = IAD10(9)
if (KTYP == 1) IADD10 = IAD10(7)
100 call dDAFILE(Lu_CIGuga,2,COP,nCOP,IADD10)
call iDAFILE(Lu_CIGuga,2,iCOP1,nCOP+1,IADD10)
LEN = ICOP1(nCOP+1)
if (LEN == 0) GO TO 100
if (LEN < 0) GO TO 200
do II=1,LEN
  IND = ICOP1(II)
  if (ICHK /= 0) GO TO 460
  if (IND /= 0) GO TO 11
  ICHK = 1
  GO TO 10
460 ICHK = 0
  ITURN = 0
  if ((IDENS == 1) .and. (IJOLD /= 0)) GO TO 20
21 ITURN = 1
  if (KTYP == 1) GO TO 9
  NK = IND
  IJOLD = NK
  NSK = NSM(NK)
  GO TO 20
9 INDI = IND
  !NI = mod(INDI,IPOW10)
  !NJ = mod(INDI/IPOW10,IPOW10)
  !NK = mod(INDI/IPOW20,IPOW10)
  NI = ibits(INDI,0,10)
  NJ = ibits(INDI,10,10)
  NK = ibits(INDI,20,10)
  NSIJ = MUL(NSM(NI),NSM(NJ))
  NSK = MUL(NSIJ,NSM(NK))
  IJ = IROW(NI)+NJ
  if (IJ == IJOLD) GO TO 20
  IJOLD = IJ
  IADR = LASTAD(NOVST+NOTT+IJ)
  do INN=1,NOB2
    FC(INN) = D0
  end do
90 call iDAFILE(Lu_TiABIJ,2,IBUFIN,LBUF2,IADR)
  LENGTH = IBUFIN(LBUF1)
  IADR = IBUFIN(LBUF2)
  if (LENGTH == 0) GO TO 91
  call SCATTER(LENGTH,FC,IBUFIN(LBUF0+1),BUFIN)
91 if (IADR /= -1) GO TO 90
  ! FORM VECTOR FK
20 NA1 = NSYS(NSK)+1
  NA2 = NSYS(NSK+1)
  INK = 0
  if (NA2 < NA1) GO TO 10
  do NA=NA1,NA2
    INK = INK+1
    NAK = IROW(LN+NA)+NK
    if (ITURN == 0) FC(NAK) = FK(INK)
    if (ITURN == 1) FK(INK) = FC(NAK)
  end do
  if (ITURN == 0) GO TO 21
  GO TO 10
11 if (INK == 0) GO TO 10
  !PAM97 ITYP = iand(IND,63)
  !PAM97 ICP2 = iand(ishft(IND,-6),8191)
  !PAM97 ICP1 = iand(ishft(IND,-19),8191)
  !ITYP = mod(IND,IPOW6)
  !ICP2 = mod(IND/IPOW6,IPOW13)
  !ICP1 = mod(IND/IPOW19,IPOW13)
  ITYP = ibits(IND,0,6)
  ICP2 = ibits(IND,6,13)
  ICP1 = ibits(IND,19,13)
  if (ITYP > 1) GO TO 12
  INDA = ICP1
  INDB = IRC(1)+ICP2
  INNY = index(INDB)+1
  if (IDENS == 1) GO TO 41
  if (INDA /= IREF0) GO TO 42
  COPI = COP(II)/sqrt(ENP(INDB))
  call DAXPY_(INK,COPI,FK,1,S(INNY),1)
  if (ITER == 1) GO TO 10
  TERM = DDOT_(INK,FK,1,C(INNY),1)
  EPP(INDB) = EPP(INDB)+COPI*TERM
  GO TO 10
42 ENPQ = (D1-THET(INDA,INDB)/D2)*(ENP(INDA)+ENP(INDB)-D1)+THET(INDA,INDB)/D2
  FACS = sqrt(ENP(INDA))*sqrt(ENP(INDB))/ENPQ
  FACW = FACS*(D2-THET(INDA,INDB))/ENPQ
  FACWA = FACW*ENP(INDA)-FACS
  FACWB = FACW*ENP(INDB)-FACS
  COPI = COP(II)*C(INDA)
  call DAXPY_(INK,COPI*FACS,FK,1,S(INNY),1)
  call DAXPY_(INK,COPI*FACWB,FK,1,W(INNY),1)
  TERM = DDOT_(INK,FK,1,C(INNY),1)
  S(INDA) = S(INDA)+COP(II)*FACS*TERM
  W(INDA) = W(INDA)+COP(II)*FACWA*TERM
  GO TO 10
41 if (INDA == IREF0) COPI = C(INDA)*COP(II)/ENP(INDB)
  ENPQ = (D1-THET(INDA,INDB)/D2)*(ENP(INDA)+ENP(INDB)-D1)+THET(INDA,INDB)/D2
  if (INDA /= IREF0) COPI = C(INDA)*COP(II)/ENPQ
  call DAXPY_(INK,COPI,C(INNY),1,FK,1)
  !write(6,654) NK,NSK,INDB
  !654 format(1X,'TYP1,NK,NSK,INDB',3I7)
  !write(6,653) (FK(I),I=1,INK)
  !653 format(1X,'FK',5F12.6)
  GO TO 10
12 if (ITER == 1) GO TO 10
  INDA = IRC(1)+ICP1
  INDB = IRC(ITYP)+ICP2
  INMY = index(INDA)+1
  INNY = index(INDB)+1
  MYSYM = JSYM(INDA)
  NYSYM = MUL(MYSYM,NSK)
  MYL = MUL(MYSYM,LSYM)
  NYL = MUL(NYSYM,LSYM)
  IFT = 0
  if (ITYP == 2) IFT = 1
  call IPO_CPF(IPOB,NVIR,MUL,NSYM,NYL,IFT)
  NVM = NVIR(MYL)
  if (IDENS == 1) GO TO 210
  ENPQ = (D1-THET(INDA,INDB)/D2)*(ENP(INDA)+ENP(INDB)-D1)+THET(INDA,INDB)/D2
  FACS = sqrt(ENP(INDA))*sqrt(ENP(INDB))/ENPQ
  FACW = FACS*(D2-THET(INDA,INDB))/ENPQ
  FACWA = FACW*ENP(INDA)-FACS
  FACWB = FACW*ENP(INDB)-FACS
  call SETZ(DBK,INK)
  call DAXPY_(INK,COP(II),FK,1,DBK,1)
  if (NYL /= 1) GO TO 25
  if (IFT == 0) call SQUAR_CPF(C(INNY+IPOB(MYL)),A,NVM)
  if (IFT == 1) call SQUARM_CPF(C(INNY+IPOB(MYL)),A,NVM)
  call SETZ(B,NVM)
  call FMMM(DBK,A,B,1,NVM,INK)
  call DAXPY_(NVM,FACS,B,1,S(INMY),1)
  call DAXPY_(NVM,FACWA,B,1,W(INMY),1)
  SIGN = D1
  if (IFT == 1) SIGN = -D1
  IOUT = INNY+IPOB(MYL)-1
  do I=1,NVM
    do J=1,I
      IOUT = IOUT+1
      TERM = DBK(I)*C(INMY+J-1)+SIGN*DBK(J)*C(INMY+I-1)
      S(IOUT) = S(IOUT)+FACS*TERM
      W(IOUT) = W(IOUT)+FACWB*TERM
    end do
    if (IFT == 1) GO TO 125
    TERM = DBK(I)*C(INMY+I-1)
    S(IOUT) = S(IOUT)-FACS*TERM
    W(IOUT) = W(IOUT)-FACWB*TERM
125 continue
  end do
  GO TO 10
25 NKM = INK*NVM
  call SETZ(B,NVM)
  if (NSK > MYL) GO TO 26
  if (IFT == 1) call VNEG_CPF(DBK,1,DBK,1,INK)
  call FMMM(DBK,C(INNY+IPOB(MYL)),B,1,NVM,INK)
  call DAXPY_(NVM,FACS,B,1,S(INMY),1)
  call DAXPY_(NVM,FACWA,B,1,W(INMY),1)
  call SETZ(B,NKM)
  call FMMM(DBK,C(INMY),B,INK,NVM,1)
  call DAXPY_(NKM,FACS,B,1,S(INNY+IPOB(MYL)),1)
  call DAXPY_(NKM,FACWB,B,1,W(INNY+IPOB(MYL)),1)
  GO TO 10
26 call FMMM(C(INNY+IPOB(NSK)),DBK,B,NVM,1,INK)
  call DAXPY_(NVM,FACS,B,1,S(INMY),1)
  call DAXPY_(NVM,FACWA,B,1,W(INMY),1)
  call SETZ(B,NKM)
  call FMMM(C(INMY),DBK,B,NVM,INK,1)
  call DAXPY_(NKM,FACS,B,1,S(INNY+IPOB(NSK)),1)
  call DAXPY_(NKM,FACWB,B,1,W(INNY+IPOB(NSK)),1)
  GO TO 10
210 call SETZ(B,INK)
  ENPQ = (D1-THET(INDA,INDB)/D2)*(ENP(INDA)+ENP(INDB)-D1)+THET(INDA,INDB)/D2
  COPI = COP(II)/ENPQ
  !write(6,652) IFT,NYL,NSK,MYL,INDA,INDB
  !652 format(1X,'TYP2',6I7)
  if (NYL /= 1) GO TO 225
  if (IFT == 0) call SQUAR_CPF(C(INNY+IPOB(MYL)),A,NVM)
  if (IFT == 1) call SQUARN_CPF(C(INNY+IPOB(MYL)),A,NVM)
  call FMMM(C(INMY),A,B,1,INK,NVM)
227 call VSMA(B,1,COPI,FK,1,FK,1,INK)
  !write(6,651) (FK(I),I=1,INK)
  !651 format(1X,'FK',5F12.6)
  GO TO 10
225 if (NSK > MYL) GO TO 226
  call FMMM(C(INNY+IPOB(MYL)),C(INMY),B,INK,1,NVM)
  if (IFT == 1) COPI = -COPI
  GO TO 227
226 call FMMM(C(INMY),C(INNY+IPOB(NSK)),B,1,INK,NVM)
  GO TO 227
10 continue
end do
GO TO 100
200 if (IDENS == 0) GO TO 201
NA1 = NSYS(NSK)+1
NA2 = NSYS(NSK+1)
INK = 0
if (NA2 < NA1) GO TO 201
do NA=NA1,NA2
  INK = INK+1
  NAK = IROW(LN+NA)+NK
  FC(NAK) = FK(INK)
end do
201 call MDSQ2(C,S,W,MUL,INDEX,JSY,NDIAG,INUM,IRC(3),LSYM,NVIRT,SQ2)
!if(IDENS == 1) write(6,876) (FC(I),I=1,NOB2)

return

end subroutine MAI
