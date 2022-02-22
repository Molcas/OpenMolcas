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

subroutine AI_CPF(JSY,INDX,C,S,FC,BUFIN,IBUFIN,A,B,FK,DBK,ENP,EPP,KTYP)
! KTYP=0  ,  (A/I)   INTEGRALS
! KTYP=1  ,  (AI/JK) INTEGRALS

use cpf_global, only: IDENS, IRC, IREF0, IROW, ITER, LASTAD, LBUF, LN, LSYM, Lu_CIGuga, Lu_TiABIJ, MUL, NDIAG, NORBT, NSM, NSYM, &
                      NSYS, NVIR, NVIRT, SQ2
use Constants, only: Zero, One
use Definitions, only: wp, iwp, r8, RtoI

implicit none
integer(kind=iwp) :: JSY(*), INDX(*), IBUFIN(*), KTYP
real(kind=wp) :: C(*), S(*), FC(*), BUFIN(*), A(*), B(*), FK(*), DBK(*), ENP(*), EPP(*)
#include "cop.fh"
integer(kind=iwp) :: I, IADR, ICHK, ICP1, ICP2, IFT, II, IJ, IJOLD, ILEN, IND, INDA, INDB, INDI, INK, INMY, INN, INNY, INUM, IOUT, &
                     IPOB(9), ITURN, ITYP, J, LBUF0, LBUF1, LBUF2, LENGTH, MYL, MYSYM, NA, NA1, NA2, NAK, NI, NJ, NK, NKM, NOB2, &
                     NOT2, NOTT, NOVST, NSIJ, NSK, NVM, NVT, NYL, NYSYM
real(kind=wp) :: COPI, SGN, TERM
logical(kind=iwp) :: Skip
integer(kind=iwp), external :: JSUNP_CPF
real(kind=r8), external :: DDOT_

NK = 0 ! dummy initialize
NSK = 0 ! dummy initialize
INUM = IRC(4)-IRC(3)
call PSQ2(C,S,MUL,INDX,JSY,NDIAG,INUM,IRC(3),LSYM,NVIRT,SQ2)
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
      if (IND == 0) then
        ICHK = 1
      else if (INK /= 0) then
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
              call DAXPY_(INK,COPI,FK,1,S(INNY),1)
              if (ITER /= 1) then
                TERM = DDOT_(INK,FK,1,C(INNY),1)
                EPP(INDB) = EPP(INDB)+COPI*TERM
              end if
            else
              COPI = COP(II)*C(INDA)
              call DAXPY_(INK,COPI,FK,1,S(INNY),1)
              TERM = DDOT_(INK,FK,1,C(INNY),1)
              S(INDA) = S(INDA)+COP(II)*TERM
            end if
          else
            if (INDA == IREF0) COPI = C(INDA)*COP(II)/ENP(INDB)
            if (INDA /= IREF0) COPI = C(INDA)*COP(II)/(sqrt(ENP(INDA))*sqrt(ENP(INDB)))
            call DAXPY_(INK,COPI,C(INNY),1,FK,1)
            !write(u6,654) NK,NSK,INDB
            !write(u6,653) (FK(I),I=1,INK)
          end if
        else if (ITER /= 1) then
          INDA = IRC(1)+ICP1
          INDB = IRC(ITYP)+ICP2
          INMY = INDX(INDA)+1
          INNY = INDX(INDB)+1
          MYSYM = JSUNP_CPF(JSY,INDA)
          NYSYM = MUL(MYSYM,NSK)
          MYL = MUL(MYSYM,LSYM)
          NYL = MUL(NYSYM,LSYM)
          IFT = 0
          if (ITYP == 2) IFT = 1
          call IPO_CPF(IPOB,NVIR,MUL,NSYM,NYL,IFT)
          NVM = NVIR(MYL)
          if (IDENS /= 1) then
            call SETZ(DBK,INK)
            call DAXPY_(INK,COP(II),FK,1,DBK,1)
            if (NYL == 1) then
              if (IFT == 0) call SQUAR_CPF(C(INNY+IPOB(MYL)),A,NVM)
              if (IFT == 1) call SQUARM_CPF(C(INNY+IPOB(MYL)),A,NVM)
              call SETZ(B,NVM)
              call FMMM(DBK,A,B,1,NVM,INK)
              call DAXPY_(NVM,One,B,1,S(INMY),1)
              SGN = One
              if (IFT == 1) SGN = -One
              IOUT = INNY+IPOB(MYL)-1
              do I=1,NVM
                do J=1,I
                  IOUT = IOUT+1
                  TERM = DBK(I)*C(INMY+J-1)+SGN*DBK(J)*C(INMY+I-1)
                  S(IOUT) = S(IOUT)+TERM
                end do
                if (IFT == 1) cycle
                TERM = DBK(I)*C(INMY+I-1)
                S(IOUT) = S(IOUT)-TERM
              end do
            else
              NKM = INK*NVM
              call SETZ(B,NVM)
              if (NSK <= MYL) then
                if (IFT == 1) call VNEG_CPF(DBK,1,DBK,1,INK)
                call FMMM(DBK,C(INNY+IPOB(MYL)),B,1,NVM,INK)
                call DAXPY_(NVM,One,B,1,S(INMY),1)
                call SETZ(B,NKM)
                call FMMM(DBK,C(INMY),B,INK,NVM,1)
                call DAXPY_(NKM,One,B,1,S(INNY+IPOB(MYL)),1)
              else
                call FMMM(C(INNY+IPOB(NSK)),DBK,B,NVM,1,INK)
                call DAXPY_(NVM,One,B,1,S(INMY),1)
                call SETZ(B,NKM)
                call FMMM(C(INMY),DBK,B,NVM,INK,1)
                call DAXPY_(NKM,One,B,1,S(INNY+IPOB(NSK)),1)
              end if
            end if
          else
            call SETZ(B,INK)
            COPI = COP(II)/(sqrt(ENP(INDA))*sqrt(ENP(INDB)))
            !write(u6,652) IFT,NYL,NSK,MYL,INDA,INDB
            if (NYL /= 1) then
              if (NSK > MYL) then
                call FMMM(C(INMY),C(INNY+IPOB(NSK)),B,1,INK,NVM)
              else
                call FMMM(C(INNY+IPOB(MYL)),C(INMY),B,INK,1,NVM)
                if (IFT == 1) COPI = -COPI
              end if
            else
              if (IFT == 0) call SQUAR_CPF(C(INNY+IPOB(MYL)),A,NVM)
              if (IFT == 1) call SQUARN_CPF(C(INNY+IPOB(MYL)),A,NVM)
              call FMMM(C(INMY),A,B,1,INK,NVM)
            end if
            call VSMA(B,1,COPI,FK,1,FK,1,INK)
            !write(u6,651) (FK(I),I=1,INK)
          end if
        end if
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
              do INN=1,NOB2
                FC(INN) = Zero
              end do
              do
                call iDAFILE(Lu_TiABIJ,2,IBUFIN,LBUF2,IADR)
                LENGTH = IBUFIN(LBUF1)
                IADR = IBUFIN(LBUF2)
                if (LENGTH /= 0) call SCATTER(LENGTH,FC,IBUFIN(LBUF0+1),BUFIN)
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
call DSQ2(C,S,MUL,INDX,JSY,NDIAG,INUM,IRC(3),LSYM,NVIRT,SQ2)

return

!651 format(1X,'FK',5F12.6)
!652 format(1X,'TYP2',6I7)
!653 format(1X,'FK',5F12.6)
!654 format(1X,'TYP1,NK,NSK,INDB',3I7)

end subroutine AI_CPF
