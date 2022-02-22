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
subroutine MIJKL(JSY,INDX,C,S,FIJKL,BUFIN,IBUFIN,W,THET,ENP,EPP,NII)

use Constants, only: Zero, One, Two, Half
use Definitions, only: wp, iwp, RtoI

implicit none
integer(kind=iwp) :: JSY(*), INDX(*), IBUFIN(*), NII
real(kind=wp) :: C(*), S(*), FIJKL(*), BUFIN(*), W(*), THET(NII,NII), ENP(*), EPP(*)
#include "cpfmcpf.fh"
#include "files_cpf.fh"
integer(kind=iwp) :: I, IADR, IC1, IC2, ICHK, IIN, ILEN, IND, INDA, INDB, INDI, INUM, IP, IVL, JP, KKBUF0, KKBUF1, KKBUF2, KP, &
                     LENGTH, LP, NA, NB, NIJ, NIJKL, NKL, NS1, NS1L
real(kind=wp) :: COPI, ENPQ, FACS, FACW, FACWA, FACWB, FINI
integer(kind=iwp), external :: JSUNP_CPF

FINI = Zero ! dummy initialize
NCONF = JSC(4)
ICHK = 0
NIJ = IROW(LN+1)
NIJKL = NIJ*(NIJ+1)/2
do I=1,NIJKL
  FIJKL(I) = Zero
end do
KKBUF0 = (RTOI*(KBUFF1+2)-2)/(RTOI+1)
KKBUF1 = RTOI*KKBUF0+KKBUF0+1
KKBUF2 = KKBUF1+1
IADR = LASTAD(1)
do
  call iDAFILE(Lu_TiABCI,2,IBUFIN,KKBUF2,IADR)
  LENGTH = IBUFIN(KKBUF1)
  IADR = IBUFIN(KKBUF2)
  if (LENGTH /= 0) call SCATTER(LENGTH,FIJKL,IBUFIN(RTOI*KKBUF0+1),BUFIN)
  if (IADR == -1) exit
end do
IADD10 = IAD10(5)
do
  call dDAFILE(Lu_CIGuga,2,COP,nCOP,IADD10)
  call iDAFILE(Lu_CIGuga,2,iCOP1,nCOP+1,IADD10)
  ILEN = ICOP1(nCOP+1)
  if (ILEN == 0) cycle
  if (ILEN < 0) exit
  do IIN=1,ILEN
    IND = ICOP1(IIN)
    if (ICHK == 0) then
      if (IND /= 0) then
        if (abs(FINI) < 1.0e-6_wp) cycle
        IVL = ibits(IND,0,6)
        IC2 = ibits(IND,6,13)
        IC1 = ibits(IND,19,13)
        COPI = COP(IIN)*FINI
        if (IVL == 0) then
          if (IC1 == IREF0) then
            COPI = COPI/sqrt(ENP(IC2))
            S(IC2) = S(IC2)+COPI
            if (ITER /= 1) EPP(IC2) = EPP(IC2)+COPI*C(IC2)
          else if (IC2 == IREF0) then
            COPI = COPI/sqrt(ENP(IC1))
            S(IC1) = S(IC1)+COPI
            if (ITER /= 1) EPP(IC1) = EPP(IC1)+COPI*C(IC1)
          else
            ENPQ = (One-THET(IC1,IC2)*Half)*(ENP(IC1)+ENP(IC2)-One)+THET(IC1,IC2)*Half
            FACS = sqrt(ENP(IC1))*sqrt(ENP(IC2))/ENPQ
            FACW = FACS*(Two-THET(IC1,IC2))/ENPQ
            FACWA = FACW*ENP(IC1)-FACS
            FACWB = FACW*ENP(IC2)-FACS
            S(IC1) = S(IC1)+FACS*COPI*C(IC2)
            S(IC2) = S(IC2)+FACS*COPI*C(IC1)
            W(IC1) = W(IC1)+FACWA*COPI*C(IC2)
            W(IC2) = W(IC2)+FACWB*COPI*C(IC1)
          end if
        else
          INDA = IRC(IVL)+IC1
          INDB = IRC(IVL)+IC2
          ENPQ = (One-THET(INDA,INDB)*Half)*(ENP(INDA)+ENP(INDB)-One)+THET(INDA,INDB)*Half
          FACS = sqrt(ENP(INDA))*sqrt(ENP(INDB))/ENPQ
          FACW = FACS*(Two-THET(INDA,INDB))/ENPQ
          FACWA = FACW*ENP(INDA)-FACS
          FACWB = FACW*ENP(INDB)-FACS
          NA = INDX(INDA)
          NB = INDX(INDB)
          NS1 = JSUNP_CPF(JSY,INDA)
          NS1L = MUL(NS1,LSYM)
          INUM = NVIR(NS1L)
          if (IVL >= 2) INUM = NNS(NS1L)
          call DAXPY_(INUM,COPI*FACS,C(NB+1),1,S(NA+1),1)
          call DAXPY_(INUM,COPI*FACS,C(NA+1),1,S(NB+1),1)
          call DAXPY_(INUM,COPI*FACWA,C(NB+1),1,W(NA+1),1)
          call DAXPY_(INUM,COPI*FACWB,C(NA+1),1,W(NB+1),1)
        end if
      else
        ICHK = 1
      end if
    else
      ICHK = 0
      INDI = IND
      IP = ibits(INDI,0,8)
      JP = ibits(INDI,8,8)
      KP = ibits(INDI,16,8)
      LP = ibits(INDI,24,8)
      NIJ = IROW(IP)+JP
      NKL = IROW(KP)+LP
      IND = NIJ*(NIJ-1)/2+NKL
      FINI = FIJKL(IND)
    end if
  end do
end do

return

end subroutine MIJKL
