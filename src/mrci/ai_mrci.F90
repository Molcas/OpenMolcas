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

!subroutine AI(INTSYM,INDX,C,S,FC,BUFIN,IBUFIN,A,B,FK,DBK,KTYP)
!subroutine AI_MRCI(INTSYM,INDX,C,S,FC,BUF,IBUF,A,B,FK,DBK,KTYP)
subroutine AI_MRCI(INTSYM,INDX,C,S,FC,A,B,FK,DBK,KTYP)

implicit real*8(A-H,O-Z)
#include "WrkSpc.fh"
#include "SysDef.fh"
#include "mrci.fh"
dimension INTSYM(*), INDX(*), C(*), S(*), FC(*), A(*), B(*), FK(*), DBK(*)
! FC(*), BUFIN(*), IBUFIN(*)
! FC(*), BUF(NBITM3), IBUF(NBITM3+2)
dimension IPOB(9)
parameter(ONE=1.0d00)
!Statement function
JSYM(L) = JSUNP(INTSYM,L)

! KTYP=0,  (A/I)   INTEGRALS
! KTYP=1,  (AI/JK) INTEGRALS

call GETMEM('BUF','ALLO','REAL',LBUF,NBITM3)
call GETMEM('IBUF','ALLO','INTE',LIBUF,NBITM3+2)

call CSCALE(INDX,INTSYM,C,SQ2)
call CSCALE(INDX,INTSYM,S,SQ2INV)
NVT = IROW(NVIRT+1)
ICHK = 0
IJOLD = 0
NK = 0
NSA = 1
NOTT = LN*(LN+1)
NOVST = LN*NVIRT+1+NVT
!PAM97 New portable code:
!PAM04 NBCMX3 = (RTOI*NBSIZ3-2)/(RTOI+1)
!PAM04 IBOFF3 = RTOI*NBCMX3
!PAM04 IBBC3 = IBOFF3+NBCMX3+1
!PAM04 IBDA3 = IBBC3+1

if (KTYP == 0) IADD10 = IAD10(9)
if (KTYP == 1) IADD10 = IAD10(7)
! READ A COP BUFFER
100 continue
call dDAFILE(LUSYMB,2,COP,nCOP,IADD10)
call iDAFILE(LUSYMB,2,iCOP1,nCOP+1,IADD10)
LEN = ICOP1(nCOP+1)
if (LEN == 0) GO TO 100
if (LEN < 0) GO TO 200
! LOOP THROUGH THE COP BUFFER:
do II=1,LEN
  IND = ICOP1(II)
  if (ICHK /= 0) then
    ! BEGIN A RATHER LONG IF-BLOCK.
    ! ICHK FLAG IS SET. THIS SIGNALS THAT PREVIOUS IND WAS 0, WHICH IS
    ! USED TO INDICATE CHANGE TO A NEW BLOCK OF COUPLING COEFFICIENTS.
    ! RESET ICHK FLAG.
    ICHK = 0
    if (KTYP == 0) then
      ! AI CASE. SAVE INTERNAL ORBITAL INDEX IN NK:
      NK = IND
      IJOLD = NK
      NSK = NSM(NK)
      NSA = NSK
      GO TO 20
    end if
    ! AIJK CASE. UNPACK INTERNAL ORBITAL INDICES INTO NI,NJ,NK:
    INDI = IND
    !NI = mod(INDI,2**10)
    !NJ = mod(INDI/2**10,2**10)
    !NK = mod(INDI/2**20,2**10)
    NI = ibits(INDI,0,10)
    NJ = ibits(INDI,10,10)
    NK = ibits(INDI,20,10)
    NSI = NSM(NI)
    NSJ = NSM(NJ)
    NSK = NSM(NK)
    NSIJ = MUL(NSI,NSJ)
    NSA = MUL(NSIJ,NSK)
    IJ = IROW(NI)+NJ
    if (IJ /= IJOLD) then
      ! NEW INTERNAL PAIR IJ. LOAD A NEW SET OF INTEGRALS INTO FC:
      IJOLD = IJ
      IADR = LASTAD(NOVST+NOTT+IJ)
      call FZERO(FC,NBTRI)

90    continue
      !PAM04 call dDAFILE(Lu_60,2,IBUFIN,NBSIZ3,IADR)
      call iDAFILE(Lu_60,2,iWORK(LIBUF),NBITM3+2,IADR)
      call dDAFILE(Lu_60,2,WORK(LBUF),NBITM3,IADR)
      LENGTH = iWORK(LIBUF+NBITM3)
      IADR = iWORK(LIBUF+NBITM3+1)
      !PAM04 LENGTH = IBUFIN(IBBC3)
      !PAM04 IADR = IBUFIN(IBDA3)
      if (LENGTH == 0) GO TO 91
      !call SCATTER(LENGTH,FC,IBUFIN(IBOFF3+1),BUFIN)
      do i=0,length-1
        !PAM04 fc(IBUFIN(IBOFF3+i)) = bufin(i)
        fc(iWORK(LIBUF+i)) = WORK(LBUF+i)
      end do
91    if (IADR /= -1) GO TO 90
    end if
20  continue
    ! FOR THIS PARTICULAR K, TRANSFER FC(NK,NA) TO ARRAY FK:
    NVIRA = NVIR(NSA)
    if (NVIRA == 0) goto 10
    do I=1,NVIRA
      NA = NVIRP(NSA)+I
      NAK = IROW(LN+NA)+NK
      FK(I) = FC(NAK)
    end do
    goto 10
  end if
  ! END OF THE LONG IF-BLOCK.
  if (IND == 0) then
    ! IND=0 SIGNALS SWITCH TO A NEW SET OF INTEGRALS.
    ICHK = 1
    GO TO 10
  end if
  ! WE ARE PROCESSING A COUPLING COEFFICIENT AS USUAL.
  if (NVIRA == 0) GO TO 10
  !ITYP = mod(IND,2**6)
  !ICP2 = mod(IND/2**6,2**13)
  !ICP1 = mod(IND/2**19,2**13)
  ITYP = ibits(IND,0,6)
  ICP2 = ibits(IND,6,13)
  ICP1 = ibits(IND,19,13)
  if (ITYP > 1) GO TO 12
  ! ITYP=1. VALENCE-SINGLES CASE.
  INDA = ICP1
  INDB = IRC(1)+ICP2
  INNY = INDX(INDB)+1
  COPI = COP(II)*C(INDA)
  call DAXPY_(NVIRA,COPI,FK,1,S(INNY),1)
  TERM = DDOT_(NVIRA,FK,1,C(INNY),1)
  S(INDA) = S(INDA)+COP(II)*TERM
  GO TO 10
12 if ((ITER == 1) .and. (IREST == 0)) GO TO 10
  INDA = IRC(1)+ICP1
  INDB = IRC(ITYP)+ICP2
  INMY = INDX(INDA)+1
  INNY = INDX(INDB)+1
  MYINTS = JSYM(INDA)
  NYINTS = MUL(MYINTS,NSA)
  MYEXTS = MUL(MYINTS,LSYM)
  NYEXTS = MUL(NYINTS,LSYM)
  IFT = 0
  if (ITYP == 2) IFT = 1
  call IPO(IPOB,NVIR,MUL,NSYM,NYEXTS,IFT)
  NVM = NVIR(MYEXTS)
  call FZERO(DBK,NVIRA)
  call DAXPY_(NVIRA,COP(II),FK,1,DBK,1)
  if (NYEXTS /= 1) GO TO 25
  if (IFT == 0) call SQUAR(C(INNY+IPOB(MYEXTS)),A,NVM)
  if (IFT == 1) call SQUARM(C(INNY+IPOB(MYEXTS)),A,NVM)
  call FZERO(B,NVM)
  call FMMM(DBK,A,B,1,NVM,NVIRA)
  call DAXPY_(NVM,ONE,B,1,S(INMY),1)
  SIGN = 1.0d00
  if (IFT == 1) SIGN = -1.0d00
  IOUT = INNY+IPOB(MYEXTS)-1
  do I=1,NVM
    do J=1,I
      IOUT = IOUT+1
      TERM = DBK(I)*C(INMY+J-1)+SIGN*DBK(J)*C(INMY+I-1)
      S(IOUT) = S(IOUT)+TERM
    end do
    if (IFT == 1) GO TO 125
    TERM = DBK(I)*C(INMY+I-1)
    S(IOUT) = S(IOUT)-TERM
125 continue
  end do
  GO TO 10
25 NKM = NVIRA*NVM
  call FZERO(B,NVM)
  if (NSA > MYEXTS) GO TO 26
  if (IFT == 1) call VNEG(DBK,1,DBK,1,NVIRA)
  call FMMM(DBK,C(INNY+IPOB(MYEXTS)),B,1,NVM,NVIRA)
  call DAXPY_(NVM,ONE,B,1,S(INMY),1)
  call FZERO(B,NKM)
  call FMMM(DBK,C(INMY),B,NVIRA,NVM,1)
  call DAXPY_(NKM,ONE,B,1,S(INNY+IPOB(MYEXTS)),1)
  GO TO 10
26 call FMMM(C(INNY+IPOB(NSA)),DBK,B,NVM,1,NVIRA)
  call DAXPY_(NVM,ONE,B,1,S(INMY),1)
  call FZERO(B,NKM)
  call FMMM(C(INMY),DBK,B,NVM,NVIRA,1)
  call DAXPY_(NKM,ONE,B,1,S(INNY+IPOB(NSA)),1)
  GO TO 10
10 continue
end do
GO TO 100
200 continue
call CSCALE(INDX,INTSYM,C,SQ2INV)
call CSCALE(INDX,INTSYM,S,SQ2)
call GETMEM('BUF','FREE','REAL',LBUF,NBITM3)
call GETMEM('IBUF','FREE','INTE',LIBUF,NBITM3+2)

return

end subroutine AI_MRCI
