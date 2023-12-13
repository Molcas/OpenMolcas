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

subroutine ABCI(INTSYM,indx,C,S,BMN,IBMN,BIAC,BICA,BUFIN)

use mrci_global, only: IADABCI, IRC, KBUFF1, LN, LSYM, Lu_70, LUSYMB, NSM, NVIRP, NVIRT, NVPAIR, SQ2, SQ2INV
use guga_util_global, only: COP, IAD10, ICOP1, nCOP
use Symmetry_Info, only: Mul
use Definitions, only: wp, iwp

#include "intent.fh"

implicit none
integer(kind=iwp), intent(in) :: INTSYM(*), indx(*)
real(kind=wp), intent(inout) :: C(*), S(*)
real(kind=wp), intent(_OUT_) :: BMN(*), BIAC(*), BICA(*), BUFIN(*)
integer(kind=iwp), intent(_OUT_) :: IBMN(*)
integer(kind=iwp) :: IAD15, IADD10, ICCB, ICHK, ICP1, ICP2, IIN, ILEN, ILOOP, INB, IND, INDA, INDB, INS, INSB, INSIN, INUMB, IOUT, &
                     IST, IT, ITYP, LB, MA, NB, NI, NSAVE, NSIB, NSLB
real(kind=wp) :: COPL, TERM
logical(kind=iwp) :: Skip
integer(kind=iwp), external :: JSUNP
real(kind=wp), external :: DDOT_

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
do
  NI = NSAVE
  IOUT = 0
  Skip = .false.
  do
    IIN = IIN+1
    if (IIN > ILEN) then
      call dDAFILE(LUSYMB,2,COP,nCOP,IADD10)
      call iDAFILE(LUSYMB,2,iCOP1,nCOP+1,IADD10)
      ILEN = ICOP1(nCOP+1)
      if (ILEN <= 0) then
        Skip = .true.
        exit
      end if
      IIN = 1
    end if
    if (ICHK /= 0) exit
    if (ICOP1(IIN) == 0) then
      ICHK = 1
    else
      IOUT = IOUT+1
      BMN(IOUT) = COP(IIN)
      IBMN(IOUT) = ICOP1(IIN)
    end if
  end do
  if (.not. Skip) then
    ICHK = 0
    NSAVE = ICOP1(IIN)
  end if
  do NB=1,NVIRT
    NSIB = MUL(NSM(LN+NB),NSM(NI))
    NSLB = MUL(NSM(LN+NB),LSYM)
    LB = NB-NVIRP(NSM(LN+NB))
    INS = NVPAIR(NSIB)
    ILOOP = 0
    do
      INSB = INS
      do
        if (INSIN >= KBUFF1) then
          call dDAFILE(Lu_70,2,BUFIN,KBUFF1,IAD15)
          INSIN = 0
        end if
        INB = KBUFF1-INSIN
        INUMB = INSB
        if (INSB > INB) INUMB = INB
        IST = INS-INSB+1
        if (ILOOP == 0) call DCOPY_(INUMB,BUFIN(INSIN+1),1,BIAC(IST),1)
        if (ILOOP == 1) call DCOPY_(INUMB,BUFIN(INSIN+1),1,BICA(IST),1)
        INSIN = INSIN+INUMB
        INSB = INSB-INUMB
        if (INSB <= 0) exit
      end do
      ILOOP = ILOOP+1
      if (ILOOP /= 1) exit
    end do
    do IT=1,IOUT
      IND = IBMN(IT)
      ICP1 = ibits(IND,19,13)
      INDA = IRC(1)+ICP1
      if (JSUNP(INTSYM,INDA) /= NSLB) cycle
      MA = indx(INDA)+LB
      ICP2 = ibits(IND,6,13)
      ITYP = ibits(IND,0,6)
      if (INS == 0) cycle
      COPL = BMN(IT)*C(MA)
      INDB = IRC(ITYP)+ICP2
      ICCB = indx(INDB)+1
      if (ITYP == 3) then
        TERM = DDOT_(INS,C(ICCB),1,BIAC,1)
        S(ICCB:ICCB+INS-1) = S(ICCB:ICCB+INS-1)+COPL*BIAC(1:INS)
      else
        TERM = DDOT_(INS,C(ICCB),1,BICA,1)
        S(ICCB:ICCB+INS-1) = S(ICCB:ICCB+INS-1)+COPL*BICA(1:INS)
      end if
      S(MA) = S(MA)+BMN(IT)*TERM
    end do
  end do
  if (ILEN < 0) exit
end do
call CSCALE(indx,INTSYM,C,SQ2INV)
call CSCALE(indx,INTSYM,S,SQ2)

return

end subroutine ABCI
