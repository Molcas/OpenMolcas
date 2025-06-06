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

subroutine MIJKL(JSY,INDX,C,S,FIJKL,BUFIN,W,THET,ENP,EPP,NII)

use, intrinsic :: iso_c_binding, only: c_f_pointer, c_loc
use cpf_global, only: IRC, IREF0, IROW, ITER, JSC, KBUFF1, LASTAD, LN, LSYM, Lu_CIGuga, Lu_TiABCI, NCONF, NNS, NVIR
use guga_util_global, only: COP, IAD10, ICOP1, nCOP
use Symmetry_Info, only: Mul
use Constants, only: Zero, One, Two, Half
use Definitions, only: wp, iwp, RtoI

#include "intent.fh"

implicit none
integer(kind=iwp), intent(in) :: JSY(*), INDX(*), NII
real(kind=wp), intent(in) :: C(*), THET(NII,NII), ENP(*)
real(kind=wp), intent(inout) :: S(*), W(*), EPP(*)
real(kind=wp), intent(_OUT_) :: FIJKL(*), BUFIN(*)
integer(kind=iwp) :: IADD10, IADR, IC1, IC2, ICHK, ILEN, IND, INDA, INDB, INDI, INUM, IP, IVL, JP, KKBUF0, KKBUF1, KKBUF2, KP, &
                     LENGTH, LP, NA, NB, NIJ, NIJKL, NKL, NS1, NS1L
real(kind=wp) :: COPI, ENPQ, FACS, FACW, FACWA, FACWB, FINI
integer(kind=iwp), external :: JSUNP

call MIJKL_INTERNAL(BUFIN)

! This is to allow type punning without an explicit interface
contains

subroutine MIJKL_INTERNAL(BUFIN)

  real(kind=wp), target, intent(_OUT_) :: BUFIN(*)
  integer(kind=iwp), pointer :: IBUFIN(:)
  integer(kind=iwp) :: IIN

  FINI = Zero ! dummy initialize
  NCONF = JSC(4)
  ICHK = 0
  NIJ = IROW(LN+1)
  NIJKL = NIJ*(NIJ+1)/2
  FIJKL(1:NIJKL) = Zero
  KKBUF0 = (RTOI*(KBUFF1+2)-2)/(RTOI+1)
  KKBUF1 = RTOI*KKBUF0+KKBUF0+1
  KKBUF2 = KKBUF1+1
  call c_f_pointer(c_loc(BUFIN),iBUFIN,[KKBUF2])
  IADR = LASTAD(1)
  do
    call iDAFILE(Lu_TiABCI,2,IBUFIN,KKBUF2,IADR)
    LENGTH = IBUFIN(KKBUF1)
    IADR = IBUFIN(KKBUF2)
    if (LENGTH /= 0) call SCATTER(LENGTH,FIJKL,IBUFIN(RTOI*KKBUF0+1:RTOI*KKBUF0+LENGTH),BUFIN)
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
            NS1 = JSUNP(JSY,INDA)
            NS1L = MUL(NS1,LSYM)
            INUM = NVIR(NS1L)
            if (IVL >= 2) INUM = NNS(NS1L)
            S(NA+1:NA+INUM) = S(NA+1:NA+INUM)+COPI*FACS*C(NB+1:NB+INUM)
            S(NB+1:NB+INUM) = S(NB+1:NB+INUM)+COPI*FACS*C(NA+1:NA+INUM)
            W(NA+1:NA+INUM) = W(NA+1:NA+INUM)+COPI*FACWA*C(NB+1:NB+INUM)
            W(NB+1:NB+INUM) = W(NB+1:NB+INUM)+COPI*FACWB*C(NA+1:NA+INUM)
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

  nullify(IBUFIN)

  return

end subroutine MIJKL_INTERNAL

end subroutine MIJKL
