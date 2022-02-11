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

subroutine ABCI_MRCI(INTSYM,indx,C,S,BMN,IBMN,BIAC,BICA,BUFIN)

use Definitions, only: wp, iwp, r8

implicit none
integer(kind=iwp) :: INTSYM(*), indx(*), IBMN(*)
real(kind=wp) :: C(*), S(*), BMN(*), BIAC(*), BICA(*), BUFIN(*)
#include "mrci.fh"
integer(kind=iwp) :: IAD15, ICCB, ICHK, ICP1, ICP2, IIN, ILEN, ILOOP, INB, IND, INDA, INDB, INS, INSB, INSIN, INUMB, IOUT, IST, &
                     IT, ITYP, LB, MA, NB, NI, NSAVE, NSIB, NSLB
real(kind=wp) :: COPL, TERM
integer(kind=iwp), external :: JSUNP
real(kind=r8), external :: DDOT_
!Statement function
integer(kind=iwp) :: JSYM, L
JSYM(L) = JSUNP(INTSYM,L)

call CSCALE(indx,INTSYM,C,SQ2)
call CSCALE(indx,INTSYM,S,SQ2INV)
ICHK = 0
INSIN = KBUFF1
IAD15 = IADABCI
IADD10 = IAD10(4)
call dDAFILE(LUSYMB,2,COP,nCOP,IADD10)
call iDAFILE(LUSYMB,2,iCOP1,nCOP+1,IADD10)
ILEN = ICOP1(nCOP+1)
IIN = 2
NSAVE = ICOP1(IIN)
100 NI = NSAVE
IOUT = 0
110 IIN = IIN+1
if (IIN <= ILEN) GO TO 15
call dDAFILE(LUSYMB,2,COP,nCOP,IADD10)
call iDAFILE(LUSYMB,2,iCOP1,nCOP+1,IADD10)
ILEN = ICOP1(nCOP+1)
if (ILEN <= 0) GO TO 5
IIN = 1
15 if (ICHK /= 0) GO TO 460
if (ICOP1(IIN) == 0) GO TO 10
IOUT = IOUT+1
BMN(IOUT) = COP(IIN)
IBMN(IOUT) = ICOP1(IIN)
GO TO 110
10 ICHK = 1
GO TO 110
460 ICHK = 0
NSAVE = ICOP1(IIN)
5 continue
do NB=1,NVIRT
  NSIB = MUL(NSM(LN+NB),NSM(NI))
  NSLB = MUL(NSM(LN+NB),LSYM)
  LB = NB-NVIRP(NSM(LN+NB))
  INS = NVPAIR(NSIB)
  ILOOP = 0
72 INSB = INS
73 if (INSIN < KBUFF1) GO TO 75
  call dDAFILE(Lu_70,2,BUFIN,KBUFF1,IAD15)
  INSIN = 0
75 INB = KBUFF1-INSIN
  INUMB = INSB
  if (INSB > INB) INUMB = INB
  IST = INS-INSB+1
  if (ILOOP == 0) call DCOPY_(INUMB,BUFIN(INSIN+1),1,BIAC(IST),1)
  if (ILOOP == 1) call DCOPY_(INUMB,BUFIN(INSIN+1),1,BICA(IST),1)
  INSIN = INSIN+INUMB
  INSB = INSB-INUMB
  if (INSB > 0) GO TO 73
  ILOOP = ILOOP+1
  if (ILOOP == 1) GO TO 72
  do IT=1,IOUT
    IND = IBMN(IT)
    !ICP1 = mod(IND/2**19,2**13)
    ICP1 = ibits(IND,19,13)
    INDA = IRC(1)+ICP1
    if (JSYM(INDA) /= NSLB) GO TO 25
    MA = indx(INDA)+LB
    !ICP2 = mod(IND/2**6,2**13)
    !ITYP = mod(IND,2**6)
    ICP2 = ibits(IND,6,13)
    ITYP = ibits(IND,0,6)
    if (INS == 0) GO TO 25
    COPL = BMN(IT)*C(MA)
    INDB = IRC(ITYP)+ICP2
    ICCB = indx(INDB)+1
    if (ITYP == 3) GO TO 26
    TERM = DDOT_(INS,C(ICCB),1,BICA,1)
    call DAXPY_(INS,COPL,BICA,1,S(ICCB),1)
    GO TO 27
26  TERM = DDOT_(INS,C(ICCB),1,BIAC,1)
    call DAXPY_(INS,COPL,BIAC,1,S(ICCB),1)
27  S(MA) = S(MA)+BMN(IT)*TERM
25  continue
  end do
end do
if (ILEN >= 0) GO TO 100
call CSCALE(indx,INTSYM,C,SQ2INV)
call CSCALE(indx,INTSYM,S,SQ2)

return

end subroutine ABCI_MRCI
