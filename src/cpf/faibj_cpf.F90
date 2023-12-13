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

subroutine FAIBJ_CPF(JSY,INDX,C,S,ABIJ,AIBJ,AJBI,BUFIN,A,B,F,FSEC,ENP,EPP)

use, intrinsic :: iso_c_binding, only: c_f_pointer, c_loc
use cpf_global, only: IRC, IREF0, IROW, ITER, LASTAD, LBUF, LN, LSYM, Lu_CIGuga, Lu_TiABIJ, NDIAG, NSM, NSYM, NVIR, NVIRT, SQ2
use guga_util_global, only: COP, IAD10, ICOP1, nCOP
use Symmetry_Info, only: Mul
use Constants, only: Zero, One, Half
use Definitions, only: wp, iwp, RtoI

#include "intent.fh"

implicit none
integer(kind=iwp), intent(in) :: JSY(*), INDX(*)
real(kind=wp), intent(inout) :: C(*), S(*), ABIJ(*), AIBJ(*), AJBI(*), FSEC(*), EPP(*)
real(kind=wp), intent(_OUT_) :: BUFIN(*), A(*), B(*), F(*)
real(kind=wp), intent(in) :: ENP(*)
integer(kind=iwp) :: IAB, IADD10, IADR, IBSYM, ICHK, ICOUP, ICOUP1, ICSYM, IFAB, IFT, IFTA, IFTB, IIN, IJ1, ILEN, ILIM, IND, INDA, &
                     INDB, INDI, INMY, INNY, INS, INUM, IPF, IPF1, IPOA(9), IPOB(9), IPOF(9), ISTAR, ITURN, ITYP, JTURN, LBUF0, &
                     LBUF1, LBUF2, LENGTH, MYL, MYSYM, NAC, NBC, NI, NJ, NOT2, NOVST, NSIJ, NVT, NYL, NYSYM
real(kind=wp) :: COPI, CPL, CPLA, CPLL, FAC, TERM
logical(kind=iwp) :: Skip
integer(kind=iwp), external :: JSUNP
real(kind=wp), external :: DDOT_

call FAIBJ_CPF_INTERNAL(BUFIN)

! This is to allow type punning without an explicit interface
contains

subroutine FAIBJ_CPF_INTERNAL(BUFIN)

  real(kind=wp), target, intent(_OUT_) :: BUFIN(*)
  integer(kind=iwp), pointer :: IBUFIN(:)
  integer(kind=iwp) :: IASYM, II

  call c_f_pointer(c_loc(BUFIN),iBUFIN,[1])

  ITYP = 0 ! dummy initialize
  ICOUP = 0 ! dummy initialize
  ICOUP1 = 0 ! dummy initialize
  INUM = IRC(4)-IRC(3)
  call PSQ2(C,S,MUL,INDX,JSY,NDIAG,INUM,IRC(3),LSYM,NVIRT,SQ2)
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
                  S(INDX(INDB)+1:INDX(INDB)+INS) = S(INDX(INDB)+1:INDX(INDB)+INS)+CPLL*FSEC(ISTAR:ISTAR+INS-1)
                  if (ITER /= 1) then
                    TERM = DDOT_(INS,C(INDX(INDB)+1),1,FSEC(ISTAR),1)
                    EPP(INDB) = EPP(INDB)+CPLL*TERM
                  end if
                else
                  COPI = CPL*C(INDA)
                  S(INDX(INDB)+1:INDX(INDB)+INS) = S(INDX(INDB)+1:INDX(INDB)+INS)+COPI*FSEC(ISTAR:ISTAR+INS-1)
                  TERM = DDOT_(INS,FSEC(ISTAR),1,C(INDX(INDB)+1),1)
                  S(INDA) = S(INDA)+CPL*TERM
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
          MYSYM = JSUNP(JSY,INDA)

          NYSYM = MUL(MYSYM,NSIJ)
          MYL = MUL(MYSYM,LSYM)
          NYL = MUL(NYSYM,LSYM)
          call IPO_CPF(IPOA,NVIR,MUL,NSYM,MYL,IFTA)
          call IPO_CPF(IPOB,NVIR,MUL,NSYM,NYL,IFTB)
          INMY = INDX(INDA)+1
          INNY = INDX(INDB)+1
          if (ITYP == 5) then
            ! DOUBLET-DOUBLET INTERACTIONS
            IIN = IPOF(MYL+1)-IPOF(MYL)
            if (IIN /= 0) then
              IPF = IPOF(MYL)
              F(1:IIN) = CPL*AIBJ(IPF+1:IPF+IIN)+CPLA*ABIJ(IPF+1:IPF+IIN)
              if (INDA == INDB) call DCOPY_(NVIR(MYL),[Zero],0,F,NVIR(MYL)+1)
              call FMMM(C(INMY),F,A,1,NVIR(NYL),NVIR(MYL))
              S(INNY:INNY+NVIR(NYL)-1) = S(INNY:INNY+NVIR(NYL)-1)+A(1:NVIR(NYL))
              if (INDA /= INDB) then
                call FMMM(F,C(INNY),A,NVIR(MYL),1,NVIR(NYL))
                S(INMY:INMY+NVIR(MYL)-1) = S(INMY:INMY+NVIR(MYL)-1)+A(1:NVIR(MYL))
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
              if (IASYM > ICSYM) then
                if (IBSYM > ICSYM) then
                  ! CASE 1 , IASYM > ICSYM AND IBSYM > ICSYM
                  IPF = IPOF(IASYM)
                  F(1:IAB) = CPL*AIBJ(IPF+1:IPF+IAB)+CPLA*ABIJ(IPF+1:IPF+IAB)
                  if (INDA == INDB) call DCOPY_(NVIR(IASYM),[Zero],0,F,NVIR(IASYM)+1)
                  call FMMM(C(INMY+IPOA(IASYM)),F,A,NVIR(ICSYM),NVIR(IBSYM),NVIR(IASYM))
                  S(INNY+IPOB(IBSYM):INNY+IPOB(IBSYM)+NBC-1) = S(INNY+IPOB(IBSYM):INNY+IPOB(IBSYM)+NBC-1)+A(1:NBC)
                  if (INDA /= INDB) then
                    IPF = IPOF(IBSYM)
                    F(1:IAB) = CPL*AJBI(IPF+1:IPF+IAB)+CPLA*ABIJ(IPF+1:IPF+IAB)
                    call FMMM(C(INNY+IPOB(IBSYM)),F,A,NVIR(ICSYM),NVIR(IASYM),NVIR(IBSYM))
                    S(INMY+IPOA(IASYM):INMY+IPOA(IASYM)+NAC-1) = S(INMY+IPOA(IASYM):INMY+IPOA(IASYM)+NAC-1)+A(1:NAC)
                  end if
                else
                  ! CASE 2 , IASYM > ICSYM AND ICSYM > OR = IBSYM
                  IPF = IPOF(IBSYM)
                  F(1:IAB) = CPL*AJBI(IPF+1:IPF+IAB)+CPLA*ABIJ(IPF+1:IPF+IAB)
                  call MTRANS(C(INMY+IPOA(IASYM)),A,NVIR(IASYM),NVIR(ICSYM))
                  call FMMM(F,A,B,NVIR(IBSYM),NVIR(ICSYM),NVIR(IASYM))
                  if (NYL == 1) then
                    A(1:NBC) = B(1:NBC)
                    if (IFTB /= 1) then
                      call SIADD(A,S(INNY+IPOB(ICSYM)),NVIR(IBSYM))
                      A(1:NBC) = Zero
                      call SQUAR(C(INNY+IPOB(IBSYM)),A,NVIR(IBSYM))
                    else
                      call TRADD(A,S(INNY+IPOB(ICSYM)),NVIR(IBSYM))
                      A(1:NBC) = Zero
                      call SQUARN(C(INNY+IPOB(IBSYM)),A,NVIR(IBSYM))
                    end if
                  else
                    if (IFTB /= 1) then
                      S(INNY+IPOB(ICSYM):INNY+IPOB(ICSYM)+NBC-1) = S(INNY+IPOB(ICSYM):INNY+IPOB(ICSYM)+NBC-1)+B(1:NBC)
                    else
                      S(INNY+IPOB(ICSYM):INNY+IPOB(ICSYM)+NBC-1) = S(INNY+IPOB(ICSYM):INNY+IPOB(ICSYM)+NBC-1)-B(1:NBC)
                    end if
                    call MTRANS(C(INNY+IPOB(ICSYM)),A,NVIR(ICSYM),NVIR(IBSYM))
                    if (IFTB == 1) A(1:NBC) = -A(1:NBC)
                  end if
                  call FMMM(A,F,B,NVIR(ICSYM),NVIR(IASYM),NVIR(IBSYM))
                  S(INMY+IPOA(IASYM):INMY+IPOA(IASYM)+NAC-1) = S(INMY+IPOA(IASYM):INMY+IPOA(IASYM)+NAC-1)+B(1:NAC)
                end if
              else
                if (IBSYM > ICSYM) then
                  ! CASE 3 , ICSYM > OR = IASYM AND IBSYM > ICSYM
                  IPF = IPOF(IASYM)
                  F(1:IAB) = CPL*AIBJ(IPF+1:IPF+IAB)+CPLA*ABIJ(IPF+1:IPF+IAB)
                  if (MYL == 1) then
                    if (IFTA == 0) call SQUAR(C(INMY+IPOA(IASYM)),A,NVIR(IASYM))
                    if (IFTA == 1) call SQUARN(C(INMY+IPOA(IASYM)),A,NVIR(IASYM))
                  else
                    call MTRANS(C(INMY+IPOA(ICSYM)),A,NVIR(ICSYM),NVIR(IASYM))
                    if (IFTA == 1) A(1:NAC) = -A(1:NAC)
                  end if
                  call FMMM(A,F,B,NVIR(ICSYM),NVIR(IBSYM),NVIR(IASYM))
                  S(INNY+IPOB(IBSYM):INNY+IPOB(IBSYM)+NBC-1) = S(INNY+IPOB(IBSYM):INNY+IPOB(IBSYM)+NBC-1)+B(1:NBC)
                  call MTRANS(C(INNY+IPOB(IBSYM)),A,NVIR(IBSYM),NVIR(ICSYM))
                  call FMMM(F,A,B,NVIR(IASYM),NVIR(ICSYM),NVIR(IBSYM))
                  if (MYL == 1) then
                    A(1:NAC) = B(1:NAC)
                    if (IFTA /= 1) then
                      call SIADD(A,S(INMY+IPOA(IASYM)),NVIR(IASYM))
                    else
                      call TRADD(A,S(INMY+IPOA(IASYM)),NVIR(IASYM))
                    end if
                    A(1:NAC) = Zero
                  else if (IFTA /= 1) then
                    S(INMY+IPOA(ICSYM):INMY+IPOA(ICSYM)+NAC-1) = S(INMY+IPOA(ICSYM):INMY+IPOA(ICSYM)+NAC-1)+B(1:NAC)
                  else
                    S(INMY+IPOA(ICSYM):INMY+IPOA(ICSYM)+NAC-1) = S(INMY+IPOA(ICSYM):INMY+IPOA(ICSYM)+NAC-1)-B(1:NAC)
                  end if
                else
                  ! CASE 4 , ICSYM > OR = IASYM AND ICSYM > OR = IBSYM
                  IPF = IPOF(IBSYM)
                  F(1:IAB) = CPL*AJBI(IPF+1:IPF+IAB)+CPLA*ABIJ(IPF+1:IPF+IAB)
                  if (INDA == INDB) call DCOPY_(NVIR(IASYM),[Zero],0,F,NVIR(IASYM)+1)
                  if (MYL == 1) then
                    if (IFTA == 0) call SQUAR(C(INMY+IPOA(IASYM)),A,NVIR(IASYM))
                    if (IFTA == 1) call SQUARM(C(INMY+IPOA(IASYM)),A,NVIR(IASYM))
                  else
                    if (IFTA == 0) call DCOPY_(NAC,C(INMY+IPOA(ICSYM)),1,A,1)
                    if (IFTA == 1) call VNEG(NAC,C(INMY+IPOA(ICSYM)),1,A,1)
                  end if
                  call FMMM(F,A,B,NVIR(IBSYM),NVIR(ICSYM),NVIR(IASYM))
                  if (NYL == 1) then
                    A(1:NBC) = B(1:NBC)
                    if (IFTB /= 1) then
                      call SIADD(A,S(INNY+IPOB(ICSYM)),NVIR(IBSYM))
                    else
                      call TRADD(A,S(INNY+IPOB(ICSYM)),NVIR(IBSYM))
                    end if
                    A(1:NBC) = Zero
                  else if (IFTB /= 1) then
                    S(INNY+IPOB(ICSYM):INNY+IPOB(ICSYM)+NBC-1) = S(INNY+IPOB(ICSYM):INNY+IPOB(ICSYM)+NBC-1)+B(1:NBC)
                  else
                    S(INNY+IPOB(ICSYM):INNY+IPOB(ICSYM)+NBC-1) = S(INNY+IPOB(ICSYM):INNY+IPOB(ICSYM)+NBC-1)-B(1:NBC)
                  end if
                  if (INDA /= INDB) then
                    IPF = IPOF(IASYM)
                    F(1:IAB) = CPL*AIBJ(IPF+1:IPF+IAB)+CPLA*ABIJ(IPF+1:IPF+IAB)
                    if (NYL == 1) then
                      if (IFTB == 0) call SQUAR(C(INNY+IPOB(IBSYM)),A,NVIR(IBSYM))
                      if (IFTB == 1) call SQUARM(C(INNY+IPOB(IBSYM)),A,NVIR(IBSYM))
                    else
                      if (IFTB == 0) call DCOPY_(NBC,C(INNY+IPOB(ICSYM)),1,A,1)
                      if (IFTB == 1) call VNEG(NBC,C(INNY+IPOB(ICSYM)),1,A,1)
                    end if
                    call FMMM(F,A,B,NVIR(IASYM),NVIR(ICSYM),NVIR(IBSYM))
                    if (MYL == 1) then
                      A(1:NAC) = B(1:NAC)
                      if (IFTA /= 1) then
                        call SIADD(A,S(INMY+IPOA(ICSYM)),NVIR(IASYM))
                      else
                        call TRADD(A,S(INMY+IPOA(ICSYM)),NVIR(IASYM))
                      end if
                      !A(1:NAC) = Zero
                    else if (IFTA /= 1) then
                      S(INMY+IPOA(ICSYM):INMY+IPOA(ICSYM)+NAC-1) = S(INMY+IPOA(ICSYM):INMY+IPOA(ICSYM)+NAC-1)+B(1:NAC)
                    else
                      S(INMY+IPOA(ICSYM):INMY+IPOA(ICSYM)+NAC-1) = S(INMY+IPOA(ICSYM):INMY+IPOA(ICSYM)+NAC-1)-B(1:NAC)
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
        ABIJ(1:ILIM) = Zero
        AIBJ(1:ILIM) = Zero
        AJBI(1:ILIM) = Zero
        if (ITER == 1) then
          Skip = .true.
        else
          ! READ (AB/IJ) INTEGRALS
          IADR = LASTAD(NOVST+IJ1)
          JTURN = 0
          Skip = .false.
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
                call SCATTER(LENGTH,ABIJ,IBUFIN(LBUF0+1:LBUF0+LENGTH),BUFIN)
              else
                call SCATTER(LENGTH,AIBJ,IBUFIN(LBUF0+1:LBUF0+LENGTH),BUFIN)
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
            if (IBSYM > IASYM) cycle
            IAB = IPOA(IASYM+1)-IPOA(IASYM)
            if (IAB == 0) cycle
            call SECORD(AIBJ(IPOF(IASYM)+1),AIBJ(IPOF(IBSYM)+1),FSEC(IIN+1),FAC,NVIR(IASYM),NVIR(IBSYM),NSIJ,IFT)
            IIN = IIN+IAB
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
                call SQUAR2(ABIJ(IPF),NVIR(IASYM))
                if (NI == NJ) call SQUAR2(AIBJ(IPF),NVIR(IASYM))
                call MTRANS(AIBJ(IPF),AJBI(IPF),NVIR(IASYM),NVIR(IBSYM))
              else
                call MTRANS(ABIJ(IPF1),ABIJ(IPF),NVIR(IASYM),NVIR(IBSYM))
                call MTRANS(AIBJ(IPF1),AJBI(IPF),NVIR(IASYM),NVIR(IBSYM))
              end if
            else
              call MTRANS(AIBJ(IPF1),AJBI(IPF),NVIR(IASYM),NVIR(IBSYM))
            end if
          end do
        end if
      end if
    end do
  end do
  call DSQ2(C,S,MUL,INDX,JSY,NDIAG,INUM,IRC(3),LSYM,NVIRT,SQ2)

  nullify(IBUFIN)

  return

end subroutine FAIBJ_CPF_INTERNAL

end subroutine FAIBJ_CPF
