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

subroutine MAI(JSY,INDX,C,S,FC,BUFIN,A,B,FK,DBK,W,THET,ENP,EPP,NII,KTYP)
! KTYP=0  ,  (A/I)   INTEGRALS
! KTYP=1  ,  (AI/JK) INTEGRALS

use, intrinsic :: iso_c_binding, only: c_f_pointer, c_loc
use cpf_global, only: IDENS, IRC, IREF0, IROW, ITER, LASTAD, LBUF, LN, LSYM, Lu_CIGuga, Lu_TiABIJ, NDIAG, NORBT, NSM, NSYM, NSYS, &
                      NVIR, NVIRT, SQ2
use guga_util_global, only: COP, IAD10, ICOP1, nCOP
use Symmetry_Info, only: Mul
use Constants, only: Zero, One, Two, Half
use Definitions, only: wp, iwp, RtoI

#include "intent.fh"

implicit none
integer(kind=iwp), intent(in) :: JSY(*), INDX(*), NII, KTYP
real(kind=wp), intent(inout) :: C(*), S(*), FK(*), W(*), EPP(*)
real(kind=wp), intent(_OUT_) :: FC(*), BUFIN(*), A(*), B(*), DBK(*)
real(kind=wp), intent(in) :: THET(NII,NII), ENP(*)
integer(kind=iwp) :: IADD10, IADR, ICHK, ICP1, ICP2, IFT, IJ, IJOLD, ILEN, IND, INDA, INDB, INDI, INK, INMY, INNY, INUM, IOUT, &
                     IPOB(9), ITURN, ITYP, LBUF0, LBUF1, LBUF2, LENGTH, MYL, MYSYM, NA1, NA2, NAK, NI, NJ, NK, NKM, NOB2, NOT2, &
                     NOTT, NOVST, NSIJ, NSK, NVM, NVT, NYL, NYSYM
real(kind=wp) :: COPI, ENPQ, FACS, FACW, FACWA, FACWB, SGN, TERM
logical(kind=iwp) :: Skip
integer(kind=iwp), external :: JSUNP
real(kind=wp), external :: DDOT_

call MAI_INTERNAL(BUFIN)

! This is to allow type punning without an explicit interface
contains

subroutine MAI_INTERNAL(BUFIN)

  real(kind=wp), target, intent(_OUT_) :: BUFIN(*)
  integer(kind=iwp), pointer :: IBUFIN(:)
  integer(kind=iwp) :: I, II, J, NA

  call c_f_pointer(c_loc(BUFIN),iBUFIN,[1])

  !if (IDENS == 1) write(u6,876) (FC(I),I=1,NOB2)
  NK = 0 ! dummy initialize
  NSK = 0 ! dummy initialize
  INUM = IRC(4)-IRC(3)
  call MPSQ2(C,S,W,MUL,INDX,JSY,NDIAG,INUM,IRC(3),LSYM,NVIRT,SQ2)
  NVT = IROW(NVIRT+1)
  ICHK = 0
  IJOLD = 0
  NOB2 = IROW(NORBT+1)
  NOT2 = IROW(LN+1)
  NOTT = 2*NOT2
  NOVST = LN*NVIRT+1+NVT
  LBUF0 = RTOI*LBUF
  LBUF1 = LBUF0+LBUF+1
  LBUF2 = LBUF1+1
  if (KTYP == 0) IADD10 = IAD10(9)
  if (KTYP == 1) IADD10 = IAD10(7)
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
          if (INK == 0) cycle
          ITYP = ibits(IND,0,6)
          ICP2 = ibits(IND,6,13)
          ICP1 = ibits(IND,19,13)
          if (ITYP <= 1) then
            INDA = ICP1
            INDB = IRC(1)+ICP2
            INNY = INDX(INDB)+1
            if (IDENS /= 1) then
              if (INDA == IREF0) then
                COPI = COP(II)/sqrt(ENP(INDB))
                S(INNY:INNY+INK-1) = S(INNY:INNY+INK-1)+COPI*FK(1:INK)
                if (ITER /= 1) then
                  TERM = DDOT_(INK,FK,1,C(INNY),1)
                  EPP(INDB) = EPP(INDB)+COPI*TERM
                end if
              else
                ENPQ = (One-THET(INDA,INDB)*Half)*(ENP(INDA)+ENP(INDB)-One)+THET(INDA,INDB)*Half
                FACS = sqrt(ENP(INDA))*sqrt(ENP(INDB))/ENPQ
                FACW = FACS*(Two-THET(INDA,INDB))/ENPQ
                FACWA = FACW*ENP(INDA)-FACS
                FACWB = FACW*ENP(INDB)-FACS
                COPI = COP(II)*C(INDA)
                S(INNY:INNY+INK-1) = S(INNY:INNY+INK-1)+COPI*FACS*FK(1:INK)
                W(INNY:INNY+INK-1) = W(INNY:INNY+INK-1)+COPI*FACWB*FK(1:INK)
                TERM = DDOT_(INK,FK,1,C(INNY),1)
                S(INDA) = S(INDA)+COP(II)*FACS*TERM
                W(INDA) = W(INDA)+COP(II)*FACWA*TERM
              end if
            else
              if (INDA == IREF0) COPI = C(INDA)*COP(II)/ENP(INDB)
              ENPQ = (One-THET(INDA,INDB)*Half)*(ENP(INDA)+ENP(INDB)-One)+THET(INDA,INDB)*Half
              if (INDA /= IREF0) COPI = C(INDA)*COP(II)/ENPQ
              FK(1:INK) = FK(1:INK)+COPI*C(INNY:INNY+INK-1)
              !write(u6,654) NK,NSK,INDB
              !write(u6,653) (FK(I),I=1,INK)
            end if
          else if (ITER /= 1) then
            INDA = IRC(1)+ICP1
            INDB = IRC(ITYP)+ICP2
            INMY = INDX(INDA)+1
            INNY = INDX(INDB)+1
            MYSYM = JSUNP(JSY,INDA)
            NYSYM = MUL(MYSYM,NSK)
            MYL = MUL(MYSYM,LSYM)
            NYL = MUL(NYSYM,LSYM)
            IFT = 0
            if (ITYP == 2) IFT = 1
            call IPO_CPF(IPOB,NVIR,MUL,NSYM,NYL,IFT)
            NVM = NVIR(MYL)
            if (IDENS /= 1) then
              ENPQ = (One-THET(INDA,INDB)*Half)*(ENP(INDA)+ENP(INDB)-One)+THET(INDA,INDB)*Half
              FACS = sqrt(ENP(INDA))*sqrt(ENP(INDB))/ENPQ
              FACW = FACS*(Two-THET(INDA,INDB))/ENPQ
              FACWA = FACW*ENP(INDA)-FACS
              FACWB = FACW*ENP(INDB)-FACS
              DBK(1:INK) = COP(II)*FK(1:INK)
              if (NYL == 1) then
                if (IFT == 0) call SQUAR(C(INNY+IPOB(MYL)),A,NVM)
                if (IFT == 1) call SQUARM(C(INNY+IPOB(MYL)),A,NVM)
                call FMMM(DBK,A,B,1,NVM,INK)
                S(INMY:INMY+NVM-1) = S(INMY:INMY+NVM-1)+FACS*B(1:NVM)
                W(INMY:INMY+NVM-1) = W(INMY:INMY+NVM-1)+FACWA*B(1:NVM)
                SGN = One
                if (IFT == 1) SGN = -One
                IOUT = INNY+IPOB(MYL)-1
                do I=1,NVM
                  do J=1,I
                    IOUT = IOUT+1
                    TERM = DBK(I)*C(INMY+J-1)+SGN*DBK(J)*C(INMY+I-1)
                    S(IOUT) = S(IOUT)+FACS*TERM
                    W(IOUT) = W(IOUT)+FACWB*TERM
                  end do
                  if (IFT /= 1) then
                    TERM = DBK(I)*C(INMY+I-1)
                    S(IOUT) = S(IOUT)-FACS*TERM
                    W(IOUT) = W(IOUT)-FACWB*TERM
                  end if
                end do
              else
                NKM = INK*NVM
                if (NSK <= MYL) then
                  if (IFT == 1) DBK(1:INK) = -DBK(1:INK)
                  I = INNY+IPOB(MYL)
                  call FMMM(DBK,C(I),B,1,NVM,INK)
                  S(INMY:INMY+NVM-1) = S(INMY:INMY+NVM-1)+FACS*B(1:NVM)
                  W(INMY:INMY+NVM-1) = W(INMY:INMY+NVM-1)+FACWA*B(1:NVM)
                  call FMMM(DBK,C(INMY),B,INK,NVM,1)
                else
                  I = INNY+IPOB(NSK)
                  call FMMM(C(I),DBK,B,NVM,1,INK)
                  S(INMY:INMY+NVM-1) = S(INMY:INMY+NVM-1)+FACS*B(1:NVM)
                  W(INMY:INMY+NVM-1) = W(INMY:INMY+NVM-1)+FACWA*B(1:NVM)
                  call FMMM(C(INMY),DBK,B,NVM,INK,1)
                end if
                S(I:I+NKM-1) = S(I:I+NKM-1)+FACS*B(1:NKM)
                W(I:I+NKM-1) = W(I:I+NKM-1)+FACWB*B(1:NKM)
              end if
            else
              ENPQ = (One-THET(INDA,INDB)*Half)*(ENP(INDA)+ENP(INDB)-One)+THET(INDA,INDB)*Half
              COPI = COP(II)/ENPQ
              !write(u6,652) IFT,NYL,NSK,MYL,INDA,INDB
              if (NYL == 1) then
                if (IFT == 0) call SQUAR(C(INNY+IPOB(MYL)),A,NVM)
                if (IFT == 1) call SQUARN(C(INNY+IPOB(MYL)),A,NVM)
                call FMMM(C(INMY),A,B,1,INK,NVM)
              else if (NSK <= MYL) then
                call FMMM(C(INNY+IPOB(MYL)),C(INMY),B,INK,1,NVM)
                if (IFT == 1) COPI = -COPI
              else
                call FMMM(C(INMY),C(INNY+IPOB(NSK)),B,1,INK,NVM)
              end if
              FK(1:INK) = FK(1:INK)+COPI*B(1:INK)
              !write(u6,651) (FK(I),I=1,INK)
            end if
          end if
        else
          ICHK = 1
        end if
      else
        ICHK = 0
        ITURN = 0
        Skip = .false.
        if ((IDENS == 1) .and. (IJOLD /= 0)) Skip = .true.
        do
          if (Skip) then
            Skip = .false.
          else
            ITURN = 1
            if (KTYP /= 1) then
              NK = IND
              IJOLD = NK
              NSK = NSM(NK)
            else
              INDI = IND
              NI = ibits(INDI,0,10)
              NJ = ibits(INDI,10,10)
              NK = ibits(INDI,20,10)
              NSIJ = MUL(NSM(NI),NSM(NJ))
              NSK = MUL(NSIJ,NSM(NK))
              IJ = IROW(NI)+NJ
              if (IJ /= IJOLD) then
                IJOLD = IJ
                IADR = LASTAD(NOVST+NOTT+IJ)
                FC(1:NOB2) = Zero
                do
                  call iDAFILE(Lu_TiABIJ,2,IBUFIN,LBUF2,IADR)
                  LENGTH = IBUFIN(LBUF1)
                  IADR = IBUFIN(LBUF2)
                  if (LENGTH /= 0) call SCATTER(LENGTH,FC,IBUFIN(LBUF0+1:LBUF0+LENGTH),BUFIN)
                  if (IADR == -1) exit
                end do
              end if
            end if
          end if
          ! FORM VECTOR FK
          NA1 = NSYS(NSK)+1
          NA2 = NSYS(NSK+1)
          INK = 0
          if (NA2 < NA1) exit
          do NA=NA1,NA2
            INK = INK+1
            NAK = IROW(LN+NA)+NK
            if (ITURN == 0) FC(NAK) = FK(INK)
            if (ITURN == 1) FK(INK) = FC(NAK)
          end do
          if (ITURN /= 0) exit
        end do
      end if
    end do
  end do
  if (IDENS /= 0) then
    NA1 = NSYS(NSK)+1
    NA2 = NSYS(NSK+1)
    INK = 0
    do NA=NA1,NA2
      INK = INK+1
      NAK = IROW(LN+NA)+NK
      FC(NAK) = FK(INK)
    end do
  end if
  call MDSQ2(C,S,W,MUL,INDX,JSY,NDIAG,INUM,IRC(3),LSYM,NVIRT,SQ2)
  !if (IDENS == 1) write(u6,876) (FC(I),I=1,NOB2)

  nullify(IBUFIN)

  return

  !651 format(1X,'FK',5F12.6)
  !652 format(1X,'TYP2',6I7)
  !653 format(1X,'FK',5F12.6)
  !654 format(1X,'TYP1,NK,NSK,INDB',3I7)
  !876 format(1X,'AI',5F12.6)

end subroutine MAI_INTERNAL

end subroutine MAI
