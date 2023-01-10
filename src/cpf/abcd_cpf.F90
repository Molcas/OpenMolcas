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

subroutine ABCD_CPF(JSY,INDX,ISAB,C,S,ACBDS,ACBDT,BUFIN)

use cpf_global, only: IPASS, IRC, IROW, JJS, KBUFF1, LN, LSYM, Lu_TiABCD, NDIAG, NSM, NSYM, NSYS, NVIRT, SQ2
use Symmetry_Info, only: Mul
use Definitions, only: wp, iwp

#include "intent.fh"

implicit none
integer(kind=iwp), intent(in) :: JSY(*), INDX(*), ISAB(*)
real(kind=wp), intent(inout) :: C(*), S(*)
real(kind=wp), intent(_OUT_) :: ACBDS(*), ACBDT(*), BUFIN(*)
integer(kind=iwp) :: IAC, IACMAX, IACMIN, IAD16, IFIN1, IFIN2, ILOOP, IN1, INB, INDA, INPS, INPT, INS, INSB, INSIN, INUM, INUMB, &
                     ISAC, IST, IST1, IST2, ISTEP, ISYM, ITAIL, NA, NC, NDMAX, NOV, NSAC, NSACL, NVT
real(kind=wp) :: TERM
real(kind=wp), external :: DDOT_

IAD16 = 0
INSIN = KBUFF1
INUM = IRC(4)-IRC(3)
call PSQ2(C,S,MUL,INDX,JSY,NDIAG,INUM,IRC(3),LSYM,NVIRT,SQ2)
NVT = IROW(NVIRT+1)
NOV = (NVT-1)/IPASS+1
IACMAX = 0
do ISTEP=1,IPASS
  IACMIN = IACMAX+1
  IACMAX = IACMAX+NOV
  if (IACMAX > NVT) IACMAX = NVT
  if (IACMIN > IACMAX) cycle
  do ISYM=1,NSYM
    IST1 = IRC(3)+JJS(ISYM+9)+1
    IFIN1 = IRC(3)+JJS(ISYM+10)
    INPS = IFIN1-IST1+1
    IST2 = IRC(2)+JJS(ISYM)+1
    IFIN2 = IRC(2)+JJS(ISYM+1)
    INPT = IFIN2-IST2+1
    ITAIL = INPS+INPT
    if (ITAIL == 0) cycle
    IN1 = -NVIRT
    do NA=1,NVIRT
      IN1 = IN1+NVIRT
      do NC=1,NA
        IAC = IROW(NA)+NC
        if (IAC < IACMIN) cycle
        if (IAC > IACMAX) cycle
        if (NA == 1) cycle
        NSAC = MUL(NSM(LN+NA),NSM(LN+NC))
        NSACL = MUL(NSAC,LSYM)
        if (NSACL /= ISYM) cycle
        ISAC = ISAB(IN1+NC)
        NDMAX = NSYS(NSM(LN+NC)+1)
        if (NDMAX > NA) NDMAX = NA
        INS = ISAB(IN1+NDMAX)
        ILOOP = 0
        do
          INSB = INS
          do
            if (INSIN >= KBUFF1) then
              call dDAFILE(Lu_TiABCD,2,BUFIN,KBUFF1,IAD16)
              INSIN = 0
            end if
            INB = KBUFF1-INSIN
            INUMB = INSB
            if (INSB > INB) INUMB = INB
            IST = INS-INSB+1
            if (ILOOP == 0) call DCOPY_(INUMB,BUFIN(INSIN+1),1,ACBDS(IST),1)
            if (ILOOP == 1) call DCOPY_(INUMB,BUFIN(INSIN+1),1,ACBDT(IST),1)
            INSIN = INSIN+INUMB
            INSB = INSB-INUMB
            if (INSB <= 0) exit
          end do
          ILOOP = ILOOP+1
          if (ILOOP /= 1) exit
        end do
        if (INPS /= 0) then
          do INDA=IST1,IFIN1
            TERM = DDOT_(INS,C(INDX(INDA)+1),1,ACBDS,1)
            S(INDX(INDA)+ISAC) = S(INDX(INDA)+ISAC)+TERM
            S(INDX(INDA)+1:INDX(INDA)+INS) = S(INDX(INDA)+1:INDX(INDA)+INS)+C(INDX(INDA)+ISAC)*ACBDS(1:INS)
          end do
        end if
        if ((INPT == 0) .or. (NA == NC)) cycle
        do INDA=IST2,IFIN2
          TERM = DDOT_(INS,C(INDX(INDA)+1),1,ACBDT,1)
          S(INDX(INDA)+ISAC) = S(INDX(INDA)+ISAC)+TERM
          S(INDX(INDA)+1:INDX(INDA)+INS) = S(INDX(INDA)+1:INDX(INDA)+INS)+C(INDX(INDA)+ISAC)*ACBDT(1:INS)
        end do
      end do
    end do
  end do
end do
call DSQ2(C,S,MUL,INDX,JSY,NDIAG,INUM,IRC(3),LSYM,NVIRT,SQ2)

return

end subroutine ABCD_CPF
