!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!***********************************************************************

subroutine FAIBJ(INTSYM,INDX,C,S,ABIJ,AIBJ,AJBI,A,B,F,FSEC)

implicit real*8(A-H,O-Z)
#include "SysDef.fh"
#include "mrci.fh"
#include "WrkSpc.fh"
dimension INTSYM(*), INDX(*), C(*), S(*), ABIJ(*), AIBJ(*), AJBI(*), A(*), B(*), F(*), FSEC(*)
dimension IPOF(9), IPOA(9), IPOB(9)
external JSUNP

call GETMEM('BUF','ALLO','REAL',LBUF,NBITM3)
call GETMEM('IBUF','ALLO','INTE',LIBUF,NBITM3+2)

!vv this code is a real compiler killer!

! POW: Unnecessary but warningstopping initializations
iTyp = -1234567
iCoup = -1234567
iCoup1 = -1234567
!call getmem('test','chec','real',ldum,ndum)

call CSCALE(INDX,INTSYM,C,SQ2)
call CSCALE(INDX,INTSYM,S,SQ2INV)
ICHK = 0
IFAB = 0
NOVST = LN*NVIRT+1+(NVIRT*(NVIRT+1))/2
NOT2 = IROW(LN+1)

IADD10 = IAD10(6)

! Long loop, reading buffers until end of buffers is signalled
! by length field holding a negative number.
300 continue
call dDAFILE(LUSYMB,2,COP,nCOP,IADD10)
call iDAFILE(LUSYMB,2,iCOP1,nCOP+1,IADD10)
LENCOP = ICOP1(nCOP+1)
if (LENCOP == 0) GO TO 300
if (LENCOP < 0) GO TO 350

! Loop over the elements of this buffer
do II=1,LENCOP
  INDCOP = ICOP1(II)
  if (ICHK /= 0) GO TO 460
  if (INDCOP /= 0) GO TO 371
  ICHK = 1
  GO TO 260

460 ICHK = 0

  ! Unpack indices NI and NJ from INDCOP
  INDI = INDCOP
  !NI = mod(INDI,2**10)
  !NJ = mod(INDI/2**10,2**10)
  NI = ibits(INDI,0,10)
  NJ = ibits(INDI,10,10)

  NSIJ = MUL(NSM(NI),NSM(NJ))
  call IPO(IPOF,NVIR,MUL,NSYM,NSIJ,-1)
  IJ1 = IROW(NI)+NJ
  ILIM = IPOF(NSYM+1)
  ! Clear matrices ABIJ, AIBJ, and AJBI.
  call FZERO(ABIJ,ILIM)
  call FZERO(AIBJ,ILIM)
  call FZERO(AJBI,ILIM)
  if ((ITER == 1) .and. (IREST == 0)) GO TO 207

  ! READ (AB/IJ) INTEGRALS

  IADR = LASTAD(NOVST+IJ1)
  JTURN = 0
201 continue

  call iDAFILE(Lu_60,2,iWORK(LIBUF),NBITM3+2,IADR)
  call dDAFILE(Lu_60,2,WORK(LBUF),NBITM3,IADR)
  LENBUF = iWORK(LIBUF+NBITM3)
  IADR = iWORK(LIBUF+NBITM3+1)
  call faibj5(LENBUF,JTURN,iWORK(LIBUF),WORK(LBUF),AIBJ,ABIJ)

  if (IADR /= -1) GO TO 201
  if (JTURN == 1) GO TO 360

  ! READ (AI/BJ) INTEGRALS

207 IADR = LASTAD(NOVST+NOT2+IJ1)
  JTURN = 1
  GO TO 201

  ! CONSTRUCT FIRST ORDER MATRICES

360 FAC = 1.0d00
  if (NI == NJ) FAC = 0.5d00
  IN = 0
  ! VV: these calls to getmem are needed to cheat some compilers.

  if (FAC < 0) call getmem('CHECK','CHEC','real',0,0)

  IFT = 0
  call faibj3(NSIJ,IFT,AIBJ,FSEC,FAC,IN,INS,IPOA,IPOF)

  if ((ITER == 1) .and. (IREST == 0)) GO TO 260
  do IASYM=1,NSYM
    NVIRA = NVIR(IASYM)
    if (NVIRA == 0) GO TO 370
    IBSYM = MUL(NSIJ,IASYM)
    NVIRB = NVIR(IBSYM)
    if (NVIRB == 0) GO TO 370
    IPF = IPOF(IASYM)+1
    IPF1 = IPOF(IBSYM)+1
    if (IASYM > IBSYM) then
      call MTRANS(AIBJ(IPF1),1,AJBI(IPF),1,NVIRA,NVIRB)
      goto 370
    end if
    if (NSIJ /= 1) then
      call MTRANS(ABIJ(IPF1),1,ABIJ(IPF),1,NVIRA,NVIRB)
      call MTRANS(AIBJ(IPF1),1,AJBI(IPF),1,NVIRA,NVIRB)
    else
      call SQUAR2(ABIJ(IPF),NVIRA)
      if (NI == NJ) call SQUAR2(AIBJ(IPF),NVIRA)
      call MTRANS(AIBJ(IPF),1,AJBI(IPF),1,NVIRA,NVIRB)
    end if
370 continue
  end do
  GO TO 260
371 continue
  if (IFAB == 1) then
    CPLA = COP(II)
    IFAB = 0
    GO TO 100
  end if
  !IFAB = mod(INDCOP,2)
  !ITURN = mod(INDCOP/2,2)
  !ITYP = mod(INDCOP/2**2,2**3)
  !ICOUP = mod(INDCOP/2**5,2**13)
  !ICOUP1 = mod(INDCOP/2**18,2**13)
  IFAB = ibits(INDCOP,0,1)
  ITURN = ibits(INDCOP,1,1)
  ITYP = ibits(INDCOP,2,3)
  ICOUP = ibits(INDCOP,5,13)
  ICOUP1 = ibits(INDCOP,18,13)
  CPL = COP(II)
  CPLA = 0.0d00
  if (IFAB /= 0) GO TO 260
  if (ITURN /= 0) goto 100
  ! FIRST ORDER INTERACTION
  INDA = ICOUP
  INDB = IRC(ITYP+1)+ICOUP1
  ISTAR = 1
  if (ITYP == 1) ISTAR = INS+1
  if (INS /= 0) then
    COPI = CPL*C(INDA)
    call DAXPY_(INS,COPI,FSEC(ISTAR),1,S(INDX(INDB)+1),1)
    TERM = DDOT_(INS,FSEC(ISTAR),1,C(INDX(INDB)+1),1)
    S(INDA) = S(INDA)+CPL*TERM
  end if
  GO TO 260

  ! INTERACTIONS BETWEEN DOUBLES AND
  ! INTERACTIONS BETWEEN SINGLES
100 if ((ITER == 1) .and. (IREST == 0)) GO TO 260

  call faibj2(IFTA,IFTB,ICOUP1,ICOUP,INDA,INDB,MYSYM,INTSYM,NYSYM,NSIJ,MYL,NYL,FACS,IPOA,IPOB,INMY,INNY,INDX,iTYP)

  if (ITYP /= 5) GO TO 71
  ! DOUBLET-DOUBLET INTERACTIONS
  IN = IPOF(MYL+1)-IPOF(MYL)
  if (IN == 0) GO TO 260
  IPF = IPOF(MYL)+1
  call DYAX(IN,CPL,AIBJ(IPF),1,F,1)
  call DAXPY_(IN,CPLA,ABIJ(IPF),1,F,1)
  if (INDA == INDB) call SETZZ(F,NVIR(MYL))
  !call DGEMTX (NVIR(MYL),NVIR(NYL),FACS,F,NVIR(MYL),C(INMY),1,S(INNY),1)
  call DGEMV_('T',NVIR(MYL),NVIR(NYL),FACS,F,NVIR(MYL),C(INMY),1,1.0d0,S(INNY),1)
  if (INDA /= INDB) then
    !call DGEMX (NVIR(MYL),NVIR(NYL),FACS,F,NVIR(MYL),C(INNY),1,S(INMY),1)
    call DGEMV_('N',NVIR(MYL),NVIR(NYL),FACS,F,NVIR(MYL),C(INNY),1,1.0d0,S(INMY),1)
  end if
  GO TO 260
  ! TRIPLET-SINGLET, SINGLET-TRIPLET,
  ! TRIPLET-TRIPLET AND SINGLET-SINGLET INTERACTIONS
71 continue

  call loop70(INTSYM,INDX,C,S,ABIJ,AIBJ,AJBI,WORK(LBUF),iWORK(LIBUF),A,B,F,FSEC,IPOF,IPOA,IPOB,MYL,NYL,INDA,INDB,INMY,INNY,IFTB, &
              IFTA,FACS,IAB,CPL,CPLA,NVIRA,NVIRC,NVIRB)

260 continue
end do
GO TO 300

350 continue
call CSCALE(INDX,INTSYM,C,SQ2INV)
call CSCALE(INDX,INTSYM,S,SQ2)
call GETMEM('BUF','FREE','REAL',LBUF,NBITM3)
call GETMEM('IBUF','FREE','INTE',LIBUF,NBITM3+2)

return

end subroutine FAIBJ
