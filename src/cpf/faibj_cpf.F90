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
subroutine FAIBJ_CPF(JSY,INDX,C,S,ABIJ,AIBJ,AJBI,BUFIN,IBUFIN,A,B,F,FSEC,ENP,EPP)

use Constants, only: Zero, One, Half
use Definitions, only: wp, iwp, r8, RtoI

implicit none
integer(kind=iwp) :: JSY(*), INDX(*), IBUFIN(*)
real(kind=wp) :: C(*), S(*), ABIJ(*), AIBJ(*), AJBI(*), BUFIN(*), A(*), B(*), F(*), FSEC(*), ENP(*), EPP(*)
#include "cpfmcpf.fh"
#include "files_cpf.fh"
integer(kind=iwp) :: IAB, IADR, IASYM, IBSYM, ICHK, ICOUP, ICOUP1, ICSYM, IFAB, IFT, IFTA, IFTB, II, IIN, IJ1, ILEN, ILIM, IND, &
                     INDA, INDB, INDI, INMY, INNY, INS, INUM, IPF, IPF1, IPOA(9), IPOB(9), IPOF(9), ISTAR, ITURN, ITYP, JTURN, &
                     LBUF0, LBUF1, LBUF2, LENGTH, MYL, MYSYM, NAC, NBC, NI, NJ, NOT2, NOVST, NSIJ, NVT, NYL, NYSYM
real(kind=wp) :: COPI, CPL, CPLA, CPLL, FAC, TERM
integer(kind=iwp), external :: JSUNP_CPF
real(kind=r8), external :: DDOT_
!parameter(IPOW5=2**5,IPOW10=2**10,IPOW18=2**18)

ITYP = 0 ! dummy initialize
ICOUP = 0 ! dummy initialize
ICOUP1 = 0 ! dummy initialize
INUM = IRC(4)-IRC(3)
call PSQ2(C,S,MUL,INDX,JSY,NDIAG,INUM,IRC(3),LSYM,NVIRT,SQ2)
NVT = IROW(NVIRT+1)
ICHK = 0
IFAB = 0
NOVST = LN*NVIRT+1+NVT
LBUF0 = RTOI*LBUF
LBUF1 = LBUF0+LBUF+1
LBUF2 = LBUF1+1
NOT2 = IROW(LN+1)
IADD10 = IAD10(6)
300 call dDAFILE(Lu_CIGuga,2,COP,nCOP,IADD10)
call iDAFILE(Lu_CIGuga,2,iCOP1,nCOP+1,IADD10)
ILEN = ICOP1(nCOP+1)
if (ILEN == 0) GO TO 300
if (ILEN < 0) GO TO 350
do II=1,ILEN
  IND = ICOP1(II)
  if (ICHK /= 0) GO TO 460
  if (IND /= 0) GO TO 371
  ICHK = 1
  GO TO 260
460 ICHK = 0
  INDI = IND
  !NI = mod(INDI,1024)
  !NJ = mod(INDI/IPOW10,1024)
  NI = ibits(INDI,0,10)
  NJ = ibits(INDI,10,10)
  NSIJ = MUL(NSM(NI),NSM(NJ))
  call IPO_CPF(IPOF,NVIR,MUL,NSYM,NSIJ,-1)
  IJ1 = IROW(NI)+NJ
  ILIM = IPOF(NSYM+1)
  call FZERO(ABIJ,ILIM)
  call FZERO(AIBJ,ILIM)
  call FZERO(AJBI,ILIM)
  if (ITER == 1) GO TO 207
  ! READ (AB/IJ) INTEGRALS
  IADR = LASTAD(NOVST+IJ1)
  JTURN = 0
201 call iDAFILE(Lu_TiABIJ,2,IBUFIN,LBUF2,IADR)
  LENGTH = IBUFIN(LBUF1)
  IADR = IBUFIN(LBUF2)
  if (LENGTH == 0) GO TO 209
  if (JTURN == 1) GO TO 203
  call SCATTER(LENGTH,ABIJ,IBUFIN(LBUF0+1),BUFIN)
  GO TO 209
203 call SCATTER(LENGTH,AIBJ,IBUFIN(LBUF0+1),BUFIN)
209 if (IADR == -1) GO TO 206
  GO TO 201
206 if (JTURN == 1) GO TO 360
  ! READ (AI/BJ) INTEGRALS
207 IADR = LASTAD(NOVST+NOT2+IJ1)
  JTURN = 1
  GO TO 201
  ! CONSTRUCT FIRST ORDER MATRICES
360 FAC = Half
  if (NI /= NJ) FAC = One
  IIN = 0
  IFT = 0
  call IPO_CPF(IPOA,NVIR,MUL,NSYM,NSIJ,IFT)
852 do IASYM=1,NSYM
    IBSYM = MUL(NSIJ,IASYM)
    if (IBSYM > IASYM) GO TO 170
    IAB = IPOA(IASYM+1)-IPOA(IASYM)
    if (IAB == 0) GO TO 170
    call SECORD(AIBJ(IPOF(IASYM)+1),AIBJ(IPOF(IBSYM)+1),FSEC(IIN+1),FAC,NVIR(IASYM),NVIR(IBSYM),NSIJ,IFT)
    IIN = IIN+IAB
170 continue
  end do
  if (IFT == 1) GO TO 853
  INS = IIN
  IFT = 1
  FAC = Zero
  GO TO 852
  ! SQUARE ABIJ
853 if (ITER == 1) GO TO 260
  do IASYM=1,NSYM
    if (NVIR(IASYM) == 0) GO TO 370
    IBSYM = MUL(NSIJ,IASYM)
    if (NVIR(IBSYM) == 0) GO TO 370
    IPF = IPOF(IASYM)+1
    IPF1 = IPOF(IBSYM)+1
    if (IASYM > IBSYM) GO TO 369
    if (NSIJ /= 1) GO TO 361
    call SQUAR2_CPF(ABIJ(IPF),NVIR(IASYM))
    if (NI /= NJ) GO TO 368
    call SQUAR2_CPF(AIBJ(IPF),NVIR(IASYM))
368 call MTRANS_CPF(AIBJ(IPF),AJBI(IPF),NVIR(IASYM),NVIR(IBSYM))
    GO TO 370
361 call MTRANS_CPF(ABIJ(IPF1),ABIJ(IPF),NVIR(IASYM),NVIR(IBSYM))
369 call MTRANS_CPF(AIBJ(IPF1),AJBI(IPF),NVIR(IASYM),NVIR(IBSYM))
370 continue
  end do
  GO TO 260
371 if (IFAB == 1) GO TO 262
  !PAM97 IFAB = iand(IND,1)
  !PAM97 ITURN = iand(ishft(IND,-1),1)
  !PAM97 ITYP = iand(ishft(IND,-2),7)
  !PAM97 ICOUP = iand(ishft(IND,-5),8191)
  !PAM97 ICOUP1 = iand(ishft(IND,-18),8191)
  !IFAB = mod(IND,2)
  !ITURN = mod(IND/2,2)
  !ITYP = mod(IND/4,8)
  !ICOUP = mod(IND/IPOW5,8192)
  !ICOUP1 = mod(IND/IPOW18,8192)
  IFAB = ibits(IND,0,1)
  ITURN = ibits(IND,1,1)
  ITYP = ibits(IND,2,3)
  ICOUP = ibits(IND,5,13)
  ICOUP1 = ibits(IND,18,13)
  CPL = COP(II)
  CPLA = Zero
  if (IFAB /= 0) GO TO 260
  if (ITURN == 0) GO TO 263
  GO TO 100
262 CPLA = COP(II)
  IFAB = 0
  GO TO 100
  ! FIRST ORDER INTERACTION
263 INDA = ICOUP
  INDB = IRC(ITYP+1)+ICOUP1
  ISTAR = 1
  if (ITYP == 1) ISTAR = INS+1
  if (INS == 0) GO TO 260
  if (INDA /= IREF0) GO TO 342
  CPLL = CPL/sqrt(ENP(INDB))
  call DAXPY_(INS,CPLL,FSEC(ISTAR),1,S(INDX(INDB)+1),1)
  if (ITER == 1) GO TO 260
  TERM = DDOT_(INS,C(INDX(INDB)+1),1,FSEC(ISTAR),1)
  EPP(INDB) = EPP(INDB)+CPLL*TERM
  GO TO 260
342 continue
  COPI = CPL*C(INDA)
  call DAXPY_(INS,COPI,FSEC(ISTAR),1,S(INDX(INDB)+1),1)
  TERM = DDOT_(INS,FSEC(ISTAR),1,C(INDX(INDB)+1),1)
  S(INDA) = S(INDA)+CPL*TERM
  GO TO 260
  ! INTERACTIONS BETWEEN DOUBLES AND
  ! INTERACTIONS BETWEEN SINGLES
100 if (ITER == 1) GO TO 260
  !call JTIME(IST)
  IFTA = 0
  IFTB = 0
  GO TO(109,110,111,112,113),ITYP
109 INDA = IRC(2)+ICOUP1
  INDB = IRC(2)+ICOUP
  IFTA = 1
  IFTB = 1
  GO TO 115
110 INDA = IRC(3)+ICOUP1
  INDB = IRC(3)+ICOUP
  GO TO 115
111 INDA = IRC(2)+ICOUP1
  INDB = IRC(3)+ICOUP
  IFTA = 1
  GO TO 115
112 INDA = IRC(3)+ICOUP1
  INDB = IRC(2)+ICOUP
  IFTB = 1
  GO TO 115
113 INDA = IRC(1)+ICOUP1
  INDB = IRC(1)+ICOUP
115 MYSYM = JSUNP_CPF(JSY,INDA)
  !JSYM(L) = JSUNP_CPF(JSY,L)

  NYSYM = MUL(MYSYM,NSIJ)
  MYL = MUL(MYSYM,LSYM)
  NYL = MUL(NYSYM,LSYM)
  call IPO_CPF(IPOA,NVIR,MUL,NSYM,MYL,IFTA)
  call IPO_CPF(IPOB,NVIR,MUL,NSYM,NYL,IFTB)
  INMY = INDX(INDA)+1
  INNY = INDX(INDB)+1
  if (ITYP /= 5) GO TO 71
  ! DOUBLET-DOUBLET INTERACTIONS
  IIN = IPOF(MYL+1)-IPOF(MYL)
  if (IIN == 0) GO TO 260
  IPF = IPOF(MYL)+1
  call SETZ(F,IIN)
  call DAXPY_(IIN,CPL,AIBJ(IPF),1,F,1)
  call DAXPY_(IIN,CPLA,ABIJ(IPF),1,F,1)
  if (INDA == INDB) call SETZZ_CPF(F,NVIR(MYL))
  call SETZ(A,NVIR(NYL))
  call FMMM(C(INMY),F,A,1,NVIR(NYL),NVIR(MYL))
  call DAXPY_(NVIR(NYL),One,A,1,S(INNY),1)
  if (INDA == INDB) GO TO 260
  call SETZ(A,NVIR(MYL))
  call FMMM(F,C(INNY),A,NVIR(MYL),1,NVIR(NYL))
  call DAXPY_(NVIR(MYL),One,A,1,S(INMY),1)
  GO TO 260
  ! TRIPLET-SINGLET , SINGLET-TRIPLET ,
  ! TRIPLET-TRIPLET AND SINGLET-SINGLET INTERACTIONS
71 do IASYM=1,NSYM
    IAB = IPOF(IASYM+1)-IPOF(IASYM)
    if (IAB == 0) GO TO 70
    ICSYM = MUL(MYL,IASYM)
    IBSYM = MUL(NYL,ICSYM)
    if ((INDA == INDB) .and. (IBSYM > IASYM)) GO TO 70
    if (NVIR(ICSYM) == 0) GO TO 70
    NAC = NVIR(IASYM)*NVIR(ICSYM)
    NBC = NVIR(IBSYM)*NVIR(ICSYM)
    if (ICSYM >= IASYM) GO TO 31
    if (ICSYM >= IBSYM) GO TO 32
    ! CASE 1 , IASYM > ICSYM AND IBSYM > ICSYM
    IPF = IPOF(IASYM)+1
    call SETZ(F,IAB)
    call DAXPY_(IAB,CPL,AIBJ(IPF),1,F,1)
    call DAXPY_(IAB,CPLA,ABIJ(IPF),1,F,1)
    if (INDA == INDB) call SETZZ_CPF(F,NVIR(IASYM))
    call SETZ(A,NBC)
    call FMMM(C(INMY+IPOA(IASYM)),F,A,NVIR(ICSYM),NVIR(IBSYM),NVIR(IASYM))
    call DAXPY_(NBC,One,A,1,S(INNY+IPOB(IBSYM)),1)
    if (INDA == INDB) GO TO 70
    IPF = IPOF(IBSYM)+1
    call SETZ(F,IAB)
    call DAXPY_(IAB,CPL,AJBI(IPF),1,F,1)
    call DAXPY_(IAB,CPLA,ABIJ(IPF),1,F,1)
    call SETZ(A,NAC)
    call FMMM(C(INNY+IPOB(IBSYM)),F,A,NVIR(ICSYM),NVIR(IASYM),NVIR(IBSYM))
    call DAXPY_(NAC,One,A,1,S(INMY+IPOA(IASYM)),1)
    GO TO 70
    ! CASE 2 , IASYM > ICSYM AND ICSYM > OR = IBSYM
32  IPF = IPOF(IBSYM)+1
    call SETZ(F,IAB)
    call DAXPY_(IAB,CPL,AJBI(IPF),1,F,1)
    call DAXPY_(IAB,CPLA,ABIJ(IPF),1,F,1)
    call MTRANS_CPF(C(INMY+IPOA(IASYM)),A,NVIR(IASYM),NVIR(ICSYM))
    call SETZ(B,NBC)
    call FMMM(F,A,B,NVIR(IBSYM),NVIR(ICSYM),NVIR(IASYM))
    if (NYL /= 1) GO TO 35
    call SETZ(A,NBC)
    call DAXPY_(NBC,One,B,1,A,1)
    if (IFTB == 1) GO TO 134
    call SIADD_CPF(A,S(INNY+IPOB(ICSYM)),NVIR(IBSYM))
    call SETZ(A,NBC)
    call SQUAR_CPF(C(INNY+IPOB(IBSYM)),A,NVIR(IBSYM))
    GO TO 36
134 call TRADD_CPF(A,S(INNY+IPOB(ICSYM)),NVIR(IBSYM))
    call SETZ(A,NBC)
    call SQUARN_CPF(C(INNY+IPOB(IBSYM)),A,NVIR(IBSYM))
    GO TO 36
35  if (IFTB == 1) GO TO 135
    call DAXPY_(NBC,One,B,1,S(INNY+IPOB(ICSYM)),1)
    GO TO 136
135 call DAXPY_(NBC,-One,B,1,S(INNY+IPOB(ICSYM)),1)
136 call MTRANS_CPF(C(INNY+IPOB(ICSYM)),A,NVIR(ICSYM),NVIR(IBSYM))
    if (IFTB == 1) call VNEG_CPF(A,1,A,1,NBC)
36  call SETZ(B,NAC)
    call FMMM(A,F,B,NVIR(ICSYM),NVIR(IASYM),NVIR(IBSYM))
    call DAXPY_(NAC,One,B,1,S(INMY+IPOA(IASYM)),1)
    GO TO 70
31  if (ICSYM >= IBSYM) GO TO 33
    ! CASE 3 , ICSYM > OR = IASYM AND IBSYM > ICSYM
    IPF = IPOF(IASYM)+1
    call SETZ(F,IAB)
    call DAXPY_(IAB,CPL,AIBJ(IPF),1,F,1)
    call DAXPY_(IAB,CPLA,ABIJ(IPF),1,F,1)
    if (MYL /= 1) GO TO 39
    if (IFTA == 0) call SQUAR_CPF(C(INMY+IPOA(IASYM)),A,NVIR(IASYM))
    if (IFTA == 1) call SQUARN_CPF(C(INMY+IPOA(IASYM)),A,NVIR(IASYM))
    GO TO 40
39  call MTRANS_CPF(C(INMY+IPOA(ICSYM)),A,NVIR(ICSYM),NVIR(IASYM))
    if (IFTA == 1) call VNEG_CPF(A,1,A,1,NAC)
40  call SETZ(B,NBC)
    call FMMM(A,F,B,NVIR(ICSYM),NVIR(IBSYM),NVIR(IASYM))
    call DAXPY_(NBC,One,B,1,S(INNY+IPOB(IBSYM)),1)
    call MTRANS_CPF(C(INNY+IPOB(IBSYM)),A,NVIR(IBSYM),NVIR(ICSYM))
    call SETZ(B,NAC)
    call FMMM(F,A,B,NVIR(IASYM),NVIR(ICSYM),NVIR(IBSYM))
    if (MYL /= 1) GO TO 46
    call SETZ(A,NAC)
    call DAXPY_(NAC,One,B,1,A,1)
    if (IFTA == 1) GO TO 146
    call SIADD_CPF(A,S(INMY+IPOA(IASYM)),NVIR(IASYM))
    call SETZ(A,NAC)
    GO TO 70
146 call TRADD_CPF(A,S(INMY+IPOA(IASYM)),NVIR(IASYM))
    call SETZ(A,NAC)
    GO TO 70
46  if (IFTA == 1) GO TO 1146
    call DAXPY_(NAC,One,B,1,S(INMY+IPOA(ICSYM)),1)
    GO TO 70
1146 call DAXPY_(NAC,-One,B,1,S(INMY+IPOA(ICSYM)),1)
    GO TO 70
    ! CASE 4 , ICSYM > OR = IASYM AND ICSYM > OR = IBSYM
33  IPF = IPOF(IBSYM)+1
    call SETZ(F,IAB)
    call DAXPY_(IAB,CPL,AJBI(IPF),1,F,1)
    call DAXPY_(IAB,CPLA,ABIJ(IPF),1,F,1)
    if (INDA == INDB) call SETZZ_CPF(F,NVIR(IASYM))
    if (MYL /= 1) GO TO 41
    if (IFTA == 0) call SQUAR_CPF(C(INMY+IPOA(IASYM)),A,NVIR(IASYM))
    if (IFTA == 1) call SQUARM_CPF(C(INMY+IPOA(IASYM)),A,NVIR(IASYM))
    GO TO 42
41  if (IFTA == 0) call DCOPY_(NAC,C(INMY+IPOA(ICSYM)),1,A,1)
    if (IFTA == 1) call VNEG_CPF(C(INMY+IPOA(ICSYM)),1,A,1,NAC)
42  call SETZ(B,NBC)
    call FMMM(F,A,B,NVIR(IBSYM),NVIR(ICSYM),NVIR(IASYM))
    if (NYL /= 1) GO TO 43
    call SETZ(A,NBC)
    call DAXPY_(NBC,One,B,1,A,1)
    if (IFTB == 1) GO TO 143
    call SIADD_CPF(A,S(INNY+IPOB(ICSYM)),NVIR(IBSYM))
    call SETZ(A,NBC)
    GO TO 44
143 call TRADD_CPF(A,S(INNY+IPOB(ICSYM)),NVIR(IBSYM))
    call SETZ(A,NBC)
    GO TO 44
43  if (IFTB == 1) GO TO 144
    call DAXPY_(NBC,One,B,1,S(INNY+IPOB(ICSYM)),1)
    GO TO 44
144 call DAXPY_(NBC,-One,B,1,S(INNY+IPOB(ICSYM)),1)
44  if (INDA == INDB) GO TO 70
    IPF = IPOF(IASYM)+1
    call SETZ(F,IAB)
    call DAXPY_(IAB,CPL,AIBJ(IPF),1,F,1)
    call DAXPY_(IAB,CPLA,ABIJ(IPF),1,F,1)
    if (NYL /= 1) GO TO 37
    if (IFTB == 0) call SQUAR_CPF(C(INNY+IPOB(IBSYM)),A,NVIR(IBSYM))
    if (IFTB == 1) call SQUARM_CPF(C(INNY+IPOB(IBSYM)),A,NVIR(IBSYM))
    GO TO 38
37  if (IFTB == 0) call DCOPY_(NBC,C(INNY+IPOB(ICSYM)),1,A,1)
    if (IFTB == 1) call VNEG_CPF(C(INNY+IPOB(ICSYM)),1,A,1,NBC)
38  call SETZ(B,NAC)
    call FMMM(F,A,B,NVIR(IASYM),NVIR(ICSYM),NVIR(IBSYM))
    if (MYL /= 1) GO TO 45
    call SETZ(A,NAC)
    call DAXPY_(NAC,One,B,1,A,1)
    if (IFTA == 1) GO TO 145
    call SIADD_CPF(A,S(INMY+IPOA(ICSYM)),NVIR(IASYM))
    !call SETZ(A,NAC)
    GO TO 70
145 call TRADD_CPF(A,S(INMY+IPOA(ICSYM)),NVIR(IASYM))
    !call SETZ(A,NAC)
    GO TO 70
45  if (IFTA == 1) GO TO 147
    call DAXPY_(NAC,One,B,1,S(INMY+IPOA(ICSYM)),1)
    GO TO 70
147 call DAXPY_(NAC,-One,B,1,S(INMY+IPOA(ICSYM)),1)
70  continue
  end do
260 continue
end do
GO TO 300
350 call DSQ2(C,S,MUL,INDX,JSY,NDIAG,INUM,IRC(3),LSYM,NVIRT,SQ2)

return

end subroutine FAIBJ_CPF
