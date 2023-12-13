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

subroutine MABCD(JSY,INDX,ISAB,C,S,ACBDS,ACBDT,BUFIN,W,THET,ENP,NII)

use cpf_global, only: IPASS, IRC, IROW, JJS, KBUFF1, LN, LSYM, Lu_TiABCD, NDIAG, NSM, NSYM, NSYS, NVIRT, SQ2
use Symmetry_Info, only: Mul
use Constants, only: One, Two, Half
use Definitions, only: wp, iwp

#include "intent.fh"

integer(kind=iwp), intent(in) :: JSY(*), INDX(*), ISAB(*), NII
real(kind=wp), intent(inout) :: C(*), S(*), W(*)
real(kind=wp), intent(_OUT_) :: ACBDS(*), ACBDT(*), BUFIN(*)
real(kind=wp), intent(in) :: THET(NII,NII), ENP(*)
integer(kind=iwp) :: I, IAC, IACMAX, IACMIN, IAD16, IFIN1, IFIN2, ILOOP, IN1, INDA, INPS, INPT, INS, INSIN, INUM, ISAC, IST1, &
                     IST2, ISTEP, ISYM, ITAIL, NA, NC, NDMAX, NOV, NSAC, NSACL, NVT
real(kind=wp) :: ENPQ, FACS, FACW, TERM
real(kind=wp), external :: DDOT_

INUM = IRC(4)-IRC(3)
call MPSQ2(C,S,W,MUL,INDX,JSY,NDIAG,INUM,IRC(3),LSYM,NVIRT,SQ2)
IAD16 = 0
INSIN = KBUFF1
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
          do I=1,INS
            if (INSIN >= KBUFF1) then
              call dDAFILE(Lu_TiABCD,2,BUFIN,KBUFF1,IAD16)
              INSIN = 0
            end if
            INSIN = INSIN+1
            if (ILOOP == 0) ACBDS(I) = BUFIN(INSIN)
            if (ILOOP == 1) ACBDT(I) = BUFIN(INSIN)
          end do
          ILOOP = ILOOP+1
          if (ILOOP /= 1) exit
        end do
        if (INPS /= 0) then
          do INDA=IST1,IFIN1
            ENPQ = (One-THET(INDA,INDA)*Half)*(ENP(INDA)+ENP(INDA)-One)+THET(INDA,INDA)*Half
            FACS = sqrt(ENP(INDA))*sqrt(ENP(INDA))/ENPQ
            FACW = (FACS*(Two-THET(INDA,INDA))/ENPQ)*ENP(INDA)-FACS
            TERM = DDOT_(INS,C(INDX(INDA)+1),1,ACBDS,1)
            S(INDX(INDA)+ISAC) = S(INDX(INDA)+ISAC)+FACS*TERM
            W(INDX(INDA)+ISAC) = W(INDX(INDA)+ISAC)+FACW*TERM
            S(INDX(INDA)+1:INDX(INDA)+INS) = S(INDX(INDA)+1:INDX(INDA)+INS)+FACS*C(INDX(INDA)+ISAC)*ACBDS(1:INS)
            W(INDX(INDA)+1:INDX(INDA)+INS) = W(INDX(INDA)+1:INDX(INDA)+INS)+FACW*C(INDX(INDA)+ISAC)*ACBDS(1:INS)
          end do
        end if
        if ((INPT == 0) .or. (NA == NC)) cycle
        do INDA=IST2,IFIN2
          ENPQ = (One-THET(INDA,INDA)*Half)*(ENP(INDA)+ENP(INDA)-One)+THET(INDA,INDA)*Half
          FACS = sqrt(ENP(INDA))*sqrt(ENP(INDA))/ENPQ
          FACW = (FACS*(Two-THET(INDA,INDA))/ENPQ)*ENP(INDA)-FACS
          TERM = DDOT_(INS,C(INDX(INDA)+1),1,ACBDT,1)
          S(INDX(INDA)+ISAC) = S(INDX(INDA)+ISAC)+FACS*TERM
          W(INDX(INDA)+ISAC) = W(INDX(INDA)+ISAC)+FACW*TERM
          S(INDX(INDA)+1:INDX(INDA)+INS) = S(INDX(INDA)+1:INDX(INDA)+INS)+FACS*C(INDX(INDA)+ISAC)*ACBDT(1:INS)
          W(INDX(INDA)+1:INDX(INDA)+INS) = W(INDX(INDA)+1:INDX(INDA)+INS)+FACW*C(INDX(INDA)+ISAC)*ACBDT(1:INS)
        end do
      end do
    end do
  end do
end do
call MDSQ2(C,S,W,MUL,INDX,JSY,NDIAG,INUM,IRC(3),LSYM,NVIRT,SQ2)

return

end subroutine MABCD
