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

subroutine MFAIBJ(JSY,INDX,C,S,ABIJ,AIBJ,AJBI,BUFIN,IBUFIN,A,B,F,FSEC,W,THET,ENP,EPP,NII)

use cpf_global, only: IRC, IREF0, IROW, ITER, LASTAD, LBUF, LN, LSYM, Lu_CIGuga, Lu_TiABIJ, MUL, NDIAG, NSM, NSYM, NVIR, NVIRT, SQ2
use Constants, only: Zero, One, Two, Half
use Definitions, only: wp, iwp, r8, RtoI

implicit none
integer(kind=iwp) :: JSY(*), INDX(*), IBUFIN(*), NII
real(kind=wp) :: C(*), S(*), ABIJ(*), AIBJ(*), AJBI(*), BUFIN(*), A(*), B(*), F(*), FSEC(*), W(*), THET(NII,NII), ENP(*), EPP(*)
#include "cop.fh"
integer(kind=iwp) :: IAB, IADR, IASYM, IBSYM, ICHK, ICOUP, ICOUP1, ICSYM, IFAB, IFT, IFTA, IFTB, II, IIN, IJ1, ILEN, ILIM, IND, &
                     INDA, INDB, INDI, INMY, INNY, INS, INUM, IPF, IPF1, IPOA(9), IPOB(9), IPOF(9), ISTAR, ITURN, ITYP, JTURN, &
                     LBUF0, LBUF1, LBUF2, LENGTH, MYL, MYSYM, NAC, NBC, NI, NJ, NOT2, NOVST, NSIJ, NVT, NYL, NYSYM
real(kind=wp) :: COPI, CPL, CPLA, CPLL, ENPQ, FAC, FACS, FACW, FACWA, FACWB, TERM
logical(kind=iwp) :: Skip
integer(kind=iwp), external :: JSUNP_CPF
real(kind=r8), external :: DDOT_

ITYP = 0 ! dummy initialize
ICOUP = 0 ! dummy initialize
ICOUP1 = 0 ! dummy initialize
INUM = IRC(4)-IRC(3)
call MPSQ2(C,S,W,MUL,INDX,JSY,NDIAG,INUM,IRC(3),LSYM,NVIRT,SQ2)
NVT = IROW(NVIRT+1)
ICHK = 0
IFAB = 0
NOVST = LN*NVIRT+1+NVT
LBUF0 = RTOI*LBUF
LBUF1 = LBUF0+LBUF+1
LBUF2 = LBUF1+1
NOT2 = IROW(LN+1)
IADD10 = IAD10(6)
do
  call dDAFILE(Lu_CIGuga,2,COP,nCOP,IADD10)
  call iDAFILE(Lu_CIGuga,2,iCOP1,nCOP+1,IADD10)
  ILEN = ICOP1(nCOP+1)
  if (ILEN == 0) cycle
  if (ILEN < 0) exit
  do II=1,ILEN
    IND = ICOP1(II)
    if (ICHK == 0) then
      if (IND /= 0) then
        if (IFAB /= 1) then
          IFAB = ibits(IND,0,1)
          ITURN = ibits(IND,1,1)
          ITYP = ibits(IND,2,3)
          ICOUP = ibits(IND,5,13)
          ICOUP1 = ibits(IND,18,13)
          CPL = COP(II)
          CPLA = Zero
          if (IFAB /= 0) cycle
          if (ITURN == 0) then
            ! FIRST ORDER INTERACTION
            INDA = ICOUP
            INDB = IRC(ITYP+1)+ICOUP1
            ISTAR = 1
            if (ITYP == 1) ISTAR = INS+1
            if (INS /= 0) then
              if (INDA == IREF0) then
                CPLL = CPL/sqrt(ENP(INDB))
                call DAXPY_(INS,CPLL,FSEC(ISTAR),1,S(INDX(INDB)+1),1)
                if (ITER /= 1) then
                  TERM = DDOT_(INS,C(INDX(INDB)+1),1,FSEC(ISTAR),1)
                  EPP(INDB) = EPP(INDB)+CPLL*TERM
                end if
              else
                ENPQ = (One-THET(INDA,INDB)*Half)*(ENP(INDA)+ENP(INDB)-One)+THET(INDA,INDB)*Half
                FACS = sqrt(ENP(INDA))*sqrt(ENP(INDB))/ENPQ
                FACW = FACS*(Two-THET(INDA,INDB))/ENPQ
                FACWA = FACW*ENP(INDA)-FACS
                FACWB = FACW*ENP(INDB)-FACS
                COPI = CPL*C(INDA)
                call DAXPY_(INS,COPI*FACS,FSEC(ISTAR),1,S(INDX(INDB)+1),1)
                call DAXPY_(INS,COPI*FACWB,FSEC(ISTAR),1,W(INDX(INDB)+1),1)
                TERM = DDOT_(INS,FSEC(ISTAR),1,C(INDX(INDB)+1),1)
                S(INDA) = S(INDA)+CPL*FACS*TERM
                W(INDA) = W(INDA)+CPL*FACWA*TERM
              end if
            end if
            cycle
          end if
        else
          CPLA = COP(II)
          IFAB = 0
        end if
        ! INTERACTIONS BETWEEN DOUBLES AND
        ! INTERACTIONS BETWEEN SINGLES
        if (ITER == 1) cycle
        !call JTIME(IST)
        IFTA = 0
        IFTB = 0
        select case (ITYP)
          case default !(1)
            INDA = IRC(2)+ICOUP1
            INDB = IRC(2)+ICOUP
            IFTA = 1
            IFTB = 1
          case (2)
            INDA = IRC(3)+ICOUP1
            INDB = IRC(3)+ICOUP
          case (3)
            INDA = IRC(2)+ICOUP1
            INDB = IRC(3)+ICOUP
            IFTA = 1
          case (4)
            INDA = IRC(3)+ICOUP1
            INDB = IRC(2)+ICOUP
            IFTB = 1
          case (5)
            INDA = IRC(1)+ICOUP1
            INDB = IRC(1)+ICOUP
        end select
        MYSYM = JSUNP_CPF(JSY,INDA)
        NYSYM = MUL(MYSYM,NSIJ)
        MYL = MUL(MYSYM,LSYM)
        NYL = MUL(NYSYM,LSYM)
        ENPQ = (One-THET(INDA,INDB)*Half)*(ENP(INDA)+ENP(INDB)-One)+THET(INDA,INDB)*Half
        FACS = sqrt(ENP(INDA))*sqrt(ENP(INDB))/ENPQ
        FACW = FACS*(Two-THET(INDA,INDB))/ENPQ
        FACWA = FACW*ENP(INDA)-FACS
        FACWB = FACW*ENP(INDB)-FACS
        call IPO_CPF(IPOA,NVIR,MUL,NSYM,MYL,IFTA)
        call IPO_CPF(IPOB,NVIR,MUL,NSYM,NYL,IFTB)
        INMY = INDX(INDA)+1
        INNY = INDX(INDB)+1
        if (ITYP == 5) then
          ! DOUBLET-DOUBLET INTERACTIONS
          IIN = IPOF(MYL+1)-IPOF(MYL)
          if (IIN /= 0) then
            IPF = IPOF(MYL)+1
            call SETZ(F,IIN)
            call DAXPY_(IIN,CPL,AIBJ(IPF),1,F,1)
            call DAXPY_(IIN,CPLA,ABIJ(IPF),1,F,1)
            if (INDA == INDB) call SETZZ_CPF(F,NVIR(MYL))
            call SETZ(A,NVIR(NYL))
            call FMMM(C(INMY),F,A,1,NVIR(NYL),NVIR(MYL))
            call DAXPY_(NVIR(NYL),FACS,A,1,S(INNY),1)
            call DAXPY_(NVIR(NYL),FACWB,A,1,W(INNY),1)
            if (INDA /= INDB) then
              call SETZ(A,NVIR(MYL))
              call FMMM(F,C(INNY),A,NVIR(MYL),1,NVIR(NYL))
              call DAXPY_(NVIR(MYL),FACS,A,1,S(INMY),1)
              call DAXPY_(NVIR(MYL),FACWA,A,1,W(INMY),1)
            end if
          end if
        else
          ! TRIPLET-SINGLET , SINGLET-TRIPLET ,
          ! TRIPLET-TRIPLET AND SINGLET-SINGLET INTERACTIONS
          do IASYM=1,NSYM
            IAB = IPOF(IASYM+1)-IPOF(IASYM)
            if (IAB == 0) cycle
            ICSYM = MUL(MYL,IASYM)
            IBSYM = MUL(NYL,ICSYM)
            if ((INDA == INDB) .and. (IBSYM > IASYM)) cycle
            if (NVIR(ICSYM) == 0) cycle
            NAC = NVIR(IASYM)*NVIR(ICSYM)
            NBC = NVIR(IBSYM)*NVIR(ICSYM)
            if (ICSYM < IASYM) then
              if (ICSYM < IBSYM) then
                ! CASE 1 , IASYM > ICSYM AND IBSYM > ICSYM
                IPF = IPOF(IASYM)+1
                call SETZ(F,IAB)
                call DAXPY_(IAB,CPL,AIBJ(IPF),1,F,1)
                call DAXPY_(IAB,CPLA,ABIJ(IPF),1,F,1)
                if (INDA == INDB) call SETZZ_CPF(F,NVIR(IASYM))
                call SETZ(A,NBC)
                call FMMM(C(INMY+IPOA(IASYM)),F,A,NVIR(ICSYM),NVIR(IBSYM),NVIR(IASYM))
                call DAXPY_(NBC,FACS,A,1,S(INNY+IPOB(IBSYM)),1)
                call DAXPY_(NBC,FACWB,A,1,W(INNY+IPOB(IBSYM)),1)
                if (INDA /= INDB) then
                  IPF = IPOF(IBSYM)+1
                  call SETZ(F,IAB)
                  call DAXPY_(IAB,CPL,AJBI(IPF),1,F,1)
                  call DAXPY_(IAB,CPLA,ABIJ(IPF),1,F,1)
                  call SETZ(A,NAC)
                  call FMMM(C(INNY+IPOB(IBSYM)),F,A,NVIR(ICSYM),NVIR(IASYM),NVIR(IBSYM))
                  call DAXPY_(NAC,FACS,A,1,S(INMY+IPOA(IASYM)),1)
                  call DAXPY_(NAC,FACWA,A,1,W(INMY+IPOA(IASYM)),1)
                end if
              else
                ! CASE 2 , IASYM > ICSYM AND ICSYM > OR = IBSYM
                IPF = IPOF(IBSYM)+1
                call SETZ(F,IAB)
                call DAXPY_(IAB,CPL,AJBI(IPF),1,F,1)
                call DAXPY_(IAB,CPLA,ABIJ(IPF),1,F,1)
                call MTRANS_CPF(C(INMY+IPOA(IASYM)),A,NVIR(IASYM),NVIR(ICSYM))
                call SETZ(B,NBC)
                call FMMM(F,A,B,NVIR(IBSYM),NVIR(ICSYM),NVIR(IASYM))
                if (NYL == 1) then
                  call SETZ(A,NBC)
                  call DAXPY_(NBC,FACS,B,1,A,1)
                  if (IFTB /= 1) then
                    call SIADD_CPF(A,S(INNY+IPOB(ICSYM)),NVIR(IBSYM))
                    call SETZ(A,NBC)
                    call DAXPY_(NBC,FACWB,B,1,A,1)
                    call SIADD_CPF(A,W(INNY+IPOB(ICSYM)),NVIR(IBSYM))
                    call SQUAR_CPF(C(INNY+IPOB(IBSYM)),A,NVIR(IBSYM))
                  else
                    call TRADD_CPF(A,S(INNY+IPOB(ICSYM)),NVIR(IBSYM))
                    call SETZ(A,NBC)
                    call DAXPY_(NBC,FACWB,B,1,A,1)
                    call TRADD_CPF(A,W(INNY+IPOB(ICSYM)),NVIR(IBSYM))
                    call SQUARN_CPF(C(INNY+IPOB(IBSYM)),A,NVIR(IBSYM))
                  end if
                else
                  if (IFTB /= 1) then
                    call DAXPY_(NBC,FACS,B,1,S(INNY+IPOB(ICSYM)),1)
                    call DAXPY_(NBC,FACWB,B,1,W(INNY+IPOB(ICSYM)),1)
                  else
                    call DAXPY_(NBC,-FACS,B,1,S(INNY+IPOB(ICSYM)),1)
                    call DAXPY_(NBC,-FACWB,B,1,W(INNY+IPOB(ICSYM)),1)
                  end if
                  call MTRANS_CPF(C(INNY+IPOB(ICSYM)),A,NVIR(ICSYM),NVIR(IBSYM))
                  if (IFTB == 1) call VNEG_CPF(A,1,A,1,NBC)
                end if
                call SETZ(B,NAC)
                call FMMM(A,F,B,NVIR(ICSYM),NVIR(IASYM),NVIR(IBSYM))
                call DAXPY_(NAC,FACS,B,1,S(INMY+IPOA(IASYM)),1)
                call DAXPY_(NAC,FACWA,B,1,W(INMY+IPOA(IASYM)),1)
              end if
            else
              if (ICSYM < IBSYM) then
                ! CASE 3 , ICSYM > OR = IASYM AND IBSYM > ICSYM
                IPF = IPOF(IASYM)+1
                call SETZ(F,IAB)
                call DAXPY_(IAB,CPL,AIBJ(IPF),1,F,1)
                call DAXPY_(IAB,CPLA,ABIJ(IPF),1,F,1)
                if (MYL == 1) then
                  if (IFTA == 0) call SQUAR_CPF(C(INMY+IPOA(IASYM)),A,NVIR(IASYM))
                  if (IFTA == 1) call SQUARN_CPF(C(INMY+IPOA(IASYM)),A,NVIR(IASYM))
                else
                  call MTRANS_CPF(C(INMY+IPOA(ICSYM)),A,NVIR(ICSYM),NVIR(IASYM))
                  if (IFTA == 1) call VNEG_CPF(A,1,A,1,NAC)
                end if
                call SETZ(B,NBC)
                call FMMM(A,F,B,NVIR(ICSYM),NVIR(IBSYM),NVIR(IASYM))
                call DAXPY_(NBC,FACS,B,1,S(INNY+IPOB(IBSYM)),1)
                call DAXPY_(NBC,FACWB,B,1,W(INNY+IPOB(IBSYM)),1)
                call MTRANS_CPF(C(INNY+IPOB(IBSYM)),A,NVIR(IBSYM),NVIR(ICSYM))
                call SETZ(B,NAC)
                call FMMM(F,A,B,NVIR(IASYM),NVIR(ICSYM),NVIR(IBSYM))
                if (MYL == 1) then
                  call SETZ(A,NAC)
                  call DAXPY_(NAC,FACS,B,1,A,1)
                  if (IFTA /= 1) then
                    call SIADD_CPF(A,S(INMY+IPOA(IASYM)),NVIR(IASYM))
                    call SETZ(A,NAC)
                    call DAXPY_(NAC,FACWA,B,1,A,1)
                    call SIADD_CPF(A,W(INMY+IPOA(IASYM)),NVIR(IASYM))
                  else
                    call TRADD_CPF(A,S(INMY+IPOA(IASYM)),NVIR(IASYM))
                    call SETZ(A,NAC)
                    call DAXPY_(NAC,FACWA,B,1,A,1)
                    call TRADD_CPF(A,W(INMY+IPOA(IASYM)),NVIR(IASYM))
                  end if
                else if (IFTA /= 1) then
                  call DAXPY_(NAC,FACS,B,1,S(INMY+IPOA(ICSYM)),1)
                  call DAXPY_(NAC,FACWA,B,1,W(INMY+IPOA(ICSYM)),1)
                else
                  call DAXPY_(NAC,-FACS,B,1,S(INMY+IPOA(ICSYM)),1)
                  call DAXPY_(NAC,-FACWA,B,1,W(INMY+IPOA(ICSYM)),1)
                end if
              else
                ! CASE 4 , ICSYM > OR = IASYM AND ICSYM > OR = IBSYM
                IPF = IPOF(IBSYM)+1
                call SETZ(F,IAB)
                call DAXPY_(IAB,CPL,AJBI(IPF),1,F,1)
                call DAXPY_(IAB,CPLA,ABIJ(IPF),1,F,1)
                if (INDA == INDB) call SETZZ_CPF(F,NVIR(IASYM))
                if (MYL == 1) then
                  if (IFTA == 0) call SQUAR_CPF(C(INMY+IPOA(IASYM)),A,NVIR(IASYM))
                  if (IFTA == 1) call SQUARM_CPF(C(INMY+IPOA(IASYM)),A,NVIR(IASYM))
                else
                  if (IFTA == 0) call DCOPY_(NAC,C(INMY+IPOA(ICSYM)),1,A,1)
                  if (IFTA == 1) call VNEG_CPF(C(INMY+IPOA(ICSYM)),1,A,1,NAC)
                end if
                call SETZ(B,NBC)
                call FMMM(F,A,B,NVIR(IBSYM),NVIR(ICSYM),NVIR(IASYM))
                if (NYL == 1) then
                  call SETZ(A,NBC)
                  call DAXPY_(NBC,FACS,B,1,A,1)
                  if (IFTB /= 1) then
                    call SIADD_CPF(A,S(INNY+IPOB(ICSYM)),NVIR(IBSYM))
                    call SETZ(A,NBC)
                    call DAXPY_(NBC,FACWB,B,1,A,1)
                    call SIADD_CPF(A,W(INNY+IPOB(ICSYM)),NVIR(IBSYM))
                  else
                    call TRADD_CPF(A,S(INNY+IPOB(ICSYM)),NVIR(IBSYM))
                    call SETZ(A,NBC)
                    call DAXPY_(NBC,FACWB,B,1,A,1)
                    call TRADD_CPF(A,W(INNY+IPOB(ICSYM)),NVIR(IBSYM))
                  end if
                else if (IFTB /= 1) then
                  call DAXPY_(NBC,FACS,B,1,S(INNY+IPOB(ICSYM)),1)
                  call DAXPY_(NBC,FACWB,B,1,W(INNY+IPOB(ICSYM)),1)
                else
                  call DAXPY_(NBC,-FACS,B,1,S(INNY+IPOB(ICSYM)),1)
                  call DAXPY_(NBC,-FACWB,B,1,W(INNY+IPOB(ICSYM)),1)
                end if
                if (INDA /= INDB) then
                  IPF = IPOF(IASYM)+1
                  call SETZ(F,IAB)
                  call DAXPY_(IAB,CPL,AIBJ(IPF),1,F,1)
                  call DAXPY_(IAB,CPLA,ABIJ(IPF),1,F,1)
                  if (NYL == 1) then
                    if (IFTB == 0) call SQUAR_CPF(C(INNY+IPOB(IBSYM)),A,NVIR(IBSYM))
                    if (IFTB == 1) call SQUARM_CPF(C(INNY+IPOB(IBSYM)),A,NVIR(IBSYM))
                  else
                    if (IFTB == 0) call DCOPY_(NBC,C(INNY+IPOB(ICSYM)),1,A,1)
                    if (IFTB == 1) call VNEG_CPF(C(INNY+IPOB(ICSYM)),1,A,1,NBC)
                  end if
                  call SETZ(B,NAC)
                  call FMMM(F,A,B,NVIR(IASYM),NVIR(ICSYM),NVIR(IBSYM))
                  if (MYL == 1) then
                    call SETZ(A,NAC)
                    call DAXPY_(NAC,FACS,B,1,A,1)
                    if (IFTA /= 1) then
                      call SIADD_CPF(A,S(INMY+IPOA(ICSYM)),NVIR(IASYM))
                      call SETZ(A,NAC)
                      call DAXPY_(NAC,FACWA,B,1,A,1)
                      call SIADD_CPF(A,W(INMY+IPOA(ICSYM)),NVIR(IASYM))
                    else
                      call TRADD_CPF(A,S(INMY+IPOA(ICSYM)),NVIR(IASYM))
                      call SETZ(A,NAC)
                      call DAXPY_(NAC,FACWA,B,1,A,1)
                      call TRADD_CPF(A,W(INMY+IPOA(ICSYM)),NVIR(IASYM))
                    end if
                  else if (IFTA /= 1) then
                    call DAXPY_(NAC,FACS,B,1,S(INMY+IPOA(ICSYM)),1)
                    call DAXPY_(NAC,FACWA,B,1,W(INMY+IPOA(ICSYM)),1)
                  else
                    call DAXPY_(NAC,-FACS,B,1,S(INMY+IPOA(ICSYM)),1)
                    call DAXPY_(NAC,-FACWA,B,1,W(INMY+IPOA(ICSYM)),1)
                  end if
                end if
              end if
            end if
          end do
        end if
      else
        ICHK = 1
      end if
    else
      ICHK = 0
      INDI = IND
      NI = ibits(INDI,0,10)
      NJ = ibits(INDI,10,10)
      NSIJ = MUL(NSM(NI),NSM(NJ))
      call IPO_CPF(IPOF,NVIR,MUL,NSYM,NSIJ,-1)
      IJ1 = IROW(NI)+NJ
      ILIM = IPOF(NSYM+1)
      call FZERO(ABIJ,ILIM)
      call FZERO(AIBJ,ILIM)
      call FZERO(AJBI,ILIM)
      if (ITER /= 1) then
        ! READ (AB/IJ) INTEGRALS
        IADR = LASTAD(NOVST+IJ1)
        JTURN = 0
        Skip = .false.
      else
        Skip = .true.
      end if
      do
        if (Skip) then
          Skip = .false.
        else
          call iDAFILE(Lu_TiABIJ,2,IBUFIN,LBUF2,IADR)
          LENGTH = IBUFIN(LBUF1)
          IADR = IBUFIN(LBUF2)
          if (LENGTH /= 0) then
            if (JTURN /= 1) then
              call SCATTER(LENGTH,ABIJ,IBUFIN(LBUF0+1),BUFIN)
            else
              call SCATTER(LENGTH,AIBJ,IBUFIN(LBUF0+1),BUFIN)
            end if
          end if
          if (IADR /= -1) cycle
          if (JTURN == 1) exit
        end if
        ! READ (AI/BJ) INTEGRALS
        IADR = LASTAD(NOVST+NOT2+IJ1)
        JTURN = 1
      end do
      ! CONSTRUCT FIRST ORDER MATRICES
      FAC = Half
      if (NI /= NJ) FAC = One
      IIN = 0
      IFT = 0
      call IPO_CPF(IPOA,NVIR,MUL,NSYM,NSIJ,IFT)
      do
        do IASYM=1,NSYM
          IBSYM = MUL(NSIJ,IASYM)
          if (IBSYM <= IASYM) then
            IAB = IPOA(IASYM+1)-IPOA(IASYM)
            if (IAB /= 0) then
              call SECORD(AIBJ(IPOF(IASYM)+1),AIBJ(IPOF(IBSYM)+1),FSEC(IIN+1),FAC,NVIR(IASYM),NVIR(IBSYM),NSIJ,IFT)
              IIN = IIN+IAB
            end if
          end if
        end do
        if (IFT == 1) exit
        INS = IIN
        IFT = 1
        FAC = Zero
      end do
      ! SQUARE ABIJ
      if (ITER /= 1) then
        do IASYM=1,NSYM
          if (NVIR(IASYM) == 0) cycle
          IBSYM = MUL(NSIJ,IASYM)
          if (NVIR(IBSYM) == 0) cycle
          IPF = IPOF(IASYM)+1
          IPF1 = IPOF(IBSYM)+1
          if (IASYM <= IBSYM) then
            if (NSIJ == 1) then
              call SQUAR2_CPF(ABIJ(IPF),NVIR(IASYM))
              if (NI == NJ) call SQUAR2_CPF(AIBJ(IPF),NVIR(IASYM))
              call MTRANS_CPF(AIBJ(IPF),AJBI(IPF),NVIR(IASYM),NVIR(IBSYM))
            else
              call MTRANS_CPF(ABIJ(IPF1),ABIJ(IPF),NVIR(IASYM),NVIR(IBSYM))
              call MTRANS_CPF(AIBJ(IPF1),AJBI(IPF),NVIR(IASYM),NVIR(IBSYM))
            end if
          else
            call MTRANS_CPF(AIBJ(IPF1),AJBI(IPF),NVIR(IASYM),NVIR(IBSYM))
          end if
        end do
      end if
    end if
  end do
end do
call MDSQ2(C,S,W,MUL,INDX,JSY,NDIAG,INUM,IRC(3),LSYM,NVIRT,SQ2)
!NCONF = JSC(4)
!write(u6,787) (S(I),I=1,NCONF)
!write(u6,786) (W(I),I=1,NCONF)

return

!786 format(1X,'W,FAIBJ',5F10.6)
!787 format(1X,'S,FAIBJ',5F10.6)

end subroutine MFAIBJ
