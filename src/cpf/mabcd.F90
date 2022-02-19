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

subroutine MABCD(JSY,INDEX,ISAB,C,S,ACBDS,ACBDT,BUFIN,W,THET,ENP,NII)

implicit real*8(A-H,O-Z)
#include "SysDef.fh"
#include "cpfmcpf.fh"
#include "files_cpf.fh"
dimension JSY(*), index(*), ISAB(*), C(*), S(*), ACBDS(*), ACBDT(*), BUFIN(*), W(*), THET(NII,NII), ENP(*)

INUM = IRC(4)-IRC(3)
call MPSQ2(C,S,W,MUL,INDEX,JSY,NDIAG,INUM,IRC(3),LSYM,NVIRT,SQ2)
IAD16 = 0
KBUFF1 = 2*9600
INSIN = KBUFF1
NVT = IROW(NVIRT+1)
NOV = (NVT-1)/IPASS+1
IACMAX = 0
do ISTEP=1,IPASS
  IACMIN = IACMAX+1
  IACMAX = IACMAX+NOV
  if (IACMAX > NVT) IACMAX = NVT
  if (IACMIN > IACMAX) GO TO 70
  do ISYM=1,NSYM
    IST1 = IRC(3)+JJS(ISYM+9)+1
    IFIN1 = IRC(3)+JJS(ISYM+10)
    INPS = IFIN1-IST1+1
    IST2 = IRC(2)+JJS(ISYM)+1
    IFIN2 = IRC(2)+JJS(ISYM+1)
    INPT = IFIN2-IST2+1
    ITAIL = INPS+INPT
    if (ITAIL == 0) GO TO 40
    IN1 = -NVIRT
    do NA=1,NVIRT
      IN1 = IN1+NVIRT
      do NC=1,NA
        IAC = IROW(NA)+NC
        if (IAC < IACMIN) GO TO 60
        if (IAC > IACMAX) GO TO 60
        if (NA == 1) GO TO 60
        NSAC = MUL(NSM(LN+NA),NSM(LN+NC))
        NSACL = MUL(NSAC,LSYM)
        if (NSACL /= ISYM) GO TO 60
        ISAC = ISAB(IN1+NC)
        NDMAX = NSYS(NSM(LN+NC)+1)
        if (NDMAX > NA) NDMAX = NA
        INS = ISAB(IN1+NDMAX)
        ILOOP = 0
72      do I=1,INS
          if (INSIN < KBUFF1) GO TO 73
          call dDAFILE(Lu_TiABCD,2,BUFIN,KBUFF1,IAD16)
          INSIN = 0
73        INSIN = INSIN+1
          if (ILOOP == 0) ACBDS(I) = BUFIN(INSIN)
          if (ILOOP == 1) ACBDT(I) = BUFIN(INSIN)
        end do
        ILOOP = ILOOP+1
        if (ILOOP == 1) GO TO 72
        if (INPS == 0) GO TO 11
        do INDA=IST1,IFIN1
          ENPQ = (D1-THET(INDA,INDA)/D2)*(ENP(INDA)+ENP(INDA)-D1)+THET(INDA,INDA)/D2
          FACS = sqrt(ENP(INDA))*sqrt(ENP(INDA))/ENPQ
          FACW = (FACS*(D2-THET(INDA,INDA))/ENPQ)*ENP(INDA)-FACS
          TERM = DDOT_(INS,C(index(INDA)+1),1,ACBDS,1)
          S(index(INDA)+ISAC) = S(index(INDA)+ISAC)+FACS*TERM
          W(index(INDA)+ISAC) = W(index(INDA)+ISAC)+FACW*TERM
          call DAXPY_(INS,FACS*C(index(INDA)+ISAC),ACBDS,1,S(index(INDA)+1),1)
          call DAXPY_(INS,FACW*C(index(INDA)+ISAC),ACBDS,1,W(index(INDA)+1),1)
        end do
11      if ((INPT == 0) .or. (NA == NC)) GO TO 60
        do INDA=IST2,IFIN2
          ENPQ = (D1-THET(INDA,INDA)/D2)*(ENP(INDA)+ENP(INDA)-D1)+THET(INDA,INDA)/D2
          FACS = sqrt(ENP(INDA))*sqrt(ENP(INDA))/ENPQ
          FACW = (FACS*(D2-THET(INDA,INDA))/ENPQ)*ENP(INDA)-FACS
          TERM = DDOT_(INS,C(index(INDA)+1),1,ACBDT,1)
          S(index(INDA)+ISAC) = S(index(INDA)+ISAC)+FACS*TERM
          W(index(INDA)+ISAC) = W(index(INDA)+ISAC)+FACW*TERM
          call DAXPY_(INS,FACS*C(index(INDA)+ISAC),ACBDT,1,S(index(INDA)+1),1)
          call DAXPY_(INS,FACW*C(index(INDA)+ISAC),ACBDT,1,W(index(INDA)+1),1)
        end do
60      continue
      end do
    end do
40  continue
  end do
70 continue
end do
call MDSQ2(C,S,W,MUL,INDEX,JSY,NDIAG,INUM,IRC(3),LSYM,NVIRT,SQ2)

return

end subroutine MABCD
