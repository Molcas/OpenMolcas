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

subroutine MAB(ICASE,JSY,INDX,C,S,FC,A,B,F,W,THET,ENP,NII)

use cpf_global, only: IDENS, IFIRST, IPRINT, IRC, IREF0, IROW, JSC, LN, LSYM, NDIAG, NORBT, NSYM, NSYS, NVIR, NVIRT, SQ2
use Symmetry_Info, only: Mul
use Constants, only: Zero, One, Two, Half
use Definitions, only: wp, iwp, u6

#include "intent.fh"

implicit none
integer(kind=iwp), intent(in) :: ICASE(*), JSY(*), INDX(*), NII
real(kind=wp), intent(inout) :: C(*), S(*), FC(*), W(*)
real(kind=wp), intent(_OUT_) :: A(*), B(*), F(*)
real(kind=wp), intent(in) :: THET(NII,NII), ENP(*)
integer(kind=iwp) :: I, IAB, IASYM, ICSYM, IFT, II1, IIA, IIC, IIN, IJ, INDA, INMY, INN, INUM, IOC(55), IPOA(9), IPF, IPOF(9), &
                     ITAIL, ITURN, JOJ, LNA, LNC, MYL, MYSYM, NA, NA1, NA2, NAA, NAB, NAC, NB, NCLIM, NOB2, NVIRA, NVIRC
real(kind=wp) :: COPI, ENPQ, FACS, FACW, RSUM, TR, TSUM
integer(kind=iwp), external :: ICUNP, JSUNP
real(kind=wp), external :: DDOT_

NAB = 0 ! dummy initialize
NOB2 = IROW(NORBT+1)
if (IPRINT >= 15) then
  write(u6,'(A,/,(10F12.6))') ' S,AB',(S(I),I=1,JSC(4))
  write(u6,'(A,/,(10F12.6))') ' W,AB',(W(I),I=1,JSC(4))
  if (IDENS == 1) write(u6,'(A,/,(10F12.6))') ' FC,AB',(FC(I),I=1,NOB2)
end if
INUM = IRC(4)-IRC(3)
call MPSQ2(C,S,W,MUL,INDX,JSY,NDIAG,INUM,IRC(3),LSYM,NVIRT,SQ2)
NCLIM = 4
if (IFIRST /= 0) NCLIM = 2
! MOVE FOCK (DENSITY) MATRIX TO F IN SYMMETRY BLOCKS
call IPO_CPF(IPOF,NVIR,MUL,NSYM,1,-1)
ITURN = 0
do
  do IASYM=1,NSYM
    IAB = IPOF(IASYM)
    NA1 = NSYS(IASYM)+1
    NA2 = NSYS(IASYM+1)
    do NA=NA1,NA2
      do NB=NA1,NA2
        IAB = IAB+1
        if (NA >= NB) NAB = IROW(LN+NA)+LN+NB
        if (NB > NA) NAB = IROW(LN+NB)+LN+NA
        if (ITURN /= 1) then
          if (IDENS == 0) F(IAB) = Zero
          if (IDENS == 1) F(IAB) = FC(NAB)
          if (NA /= NB) F(IAB) = FC(NAB)
        else
          if (NA < NB) FC(NAB) = F(IAB)
        end if
      end do
    end do
  end do
  if (ITURN /= 0) then
    TR = Zero
    IJ = 0
    do I=1,NORBT
      IJ = IJ+I
      TR = TR+FC(IJ)
    end do
    if (iPrint >= 15) write(u6,310) TR
    exit
  end if
  II1 = 0
  ITAIL = IRC(NCLIM)
  do INDA=1,ITAIL
    if (IDENS /= 0) then
      do I=1,LN
        II1 = II1+1
        JOJ = ICUNP(ICASE,II1)
        if (JOJ > 1) JOJ = JOJ-1
        IOC(I) = JOJ
      end do
    end if
    if (INDA <= IRC(1)) then
      if ((IDENS == 0) .or. (INDA == IREF0)) cycle
      ENPQ = (One-THET(INDA,INDA)*Half)*(ENP(INDA)+ENP(INDA)-One)+THET(INDA,INDA)*Half
      TSUM = C(INDA)*C(INDA)/ENPQ
    else
      MYSYM = JSUNP(JSY,INDA)
      MYL = MUL(MYSYM,LSYM)
      INMY = INDX(INDA)+1
      ENPQ = (One-THET(INDA,INDA)*Half)*(ENP(INDA)+ENP(INDA)-One)+THET(INDA,INDA)*Half
      FACS = sqrt(ENP(INDA))*sqrt(ENP(INDA))/ENPQ
      FACW = (FACS*(Two-THET(INDA,INDA))/ENPQ)*ENP(INDA)-FACS
      if (INDA > IRC(2)) then
        ! TRIPLET-TRIPLET AND SINGLET-SINGLET INTERACTIONS
        IFT = 1
        if (INDA > IRC(3)) IFT = 0
        call IPO_CPF(IPOA,NVIR,MUL,NSYM,MYL,IFT)
        IIN = 0
        TSUM = Zero
        do IASYM=1,NSYM
          IAB = IPOF(IASYM+1)-IPOF(IASYM)
          if (IAB == 0) cycle
          ICSYM = MUL(MYL,IASYM)
          if (NVIR(ICSYM) == 0) cycle
          if (IDENS /= 1) then
            if (MYL == 1) then
              if (IFT == 0) call SQUAR(C(INMY+IPOA(IASYM)),A,NVIR(IASYM))
              !if (IFT == 1) call SQUARN(C(INMY+IPOA(IASYM)),A,NVIR(IASYM))
              if (IFT == 1) call SQUARM(C(INMY+IPOA(IASYM)),A,NVIR(IASYM))
              NAA = NVIR(IASYM)*NVIR(IASYM)
              call FMMM(F(IPOF(IASYM)+1),A,B,NVIR(IASYM),NVIR(IASYM),NVIR(IASYM))
              A(1:NAA) = FACS*B(1:NAA)
              if (IFT /= 1) then
                call SIADD(A,S(INMY+IPOA(IASYM)),NVIR(IASYM))
                A(1:NAA) = FACW*B(1:NAA)
                call SIADD(A,W(INMY+IPOA(IASYM)),NVIR(IASYM))
              else
                call TRADD(A,S(INMY+IPOA(IASYM)),NVIR(IASYM))
                A(1:NAA) = FACW*B(1:NAA)
                call TRADD(A,W(INMY+IPOA(IASYM)),NVIR(IASYM))
              end if
            else
              NAC = NVIR(IASYM)*NVIR(ICSYM)
              if (IASYM <= ICSYM) then
                I = INMY+IPOA(ICSYM)
                call FMMM(F(IPOF(IASYM)+1),C(I),A,NVIR(IASYM),NVIR(ICSYM),NVIR(IASYM))
              else
                I = INMY+IPOA(IASYM)
                call FMMM(C(I),F(IPOF(IASYM)+1),A,NVIR(ICSYM),NVIR(IASYM),NVIR(IASYM))
              end if
              S(I:I+NAC-1) = S(I:I+NAC-1)+FACS*A(1:NAC)
              W(I:I+NAC-1) = W(I:I+NAC-1)+FACW*A(1:NAC)
            end if
          else
            if (MYL == 1) then
              if (IFT == 0) call SQUAR(C(INMY+IPOA(IASYM)),A,NVIR(IASYM))
              if (IFT == 1) call SQUARM(C(INMY+IPOA(IASYM)),A,NVIR(IASYM))
            else if (IASYM <= ICSYM) then
              NAC = NVIR(IASYM)*NVIR(ICSYM)
              if (IFT == 0) call DCOPY_(NAC,C(INMY+IPOA(ICSYM)),1,A,1)
              if (IFT == 1) call VNEG(NAC,C(INMY+IPOA(ICSYM)),1,A,1)
            else
              call MTRANS(C(INMY+IPOA(IASYM)),A,NVIR(IASYM),NVIR(ICSYM))
            end if
            call FMUL2(A,A,B,NVIR(IASYM),NVIR(IASYM),NVIR(ICSYM))
            IPF = IPOF(IASYM)+1
            ENPQ = (One-THET(INDA,INDA)*Half)*(ENP(INDA)+ENP(INDA)-One)+THET(INDA,INDA)*Half
            COPI = One/ENPQ
            F(IPF:IPF+IAB-1) = F(IPF:IPF+IAB-1)+COPI*B(1:IAB)
            NVIRA = NVIR(IASYM)
            NVIRC = NVIR(ICSYM)
            INN = 1
            LNC = LN+NSYS(ICSYM)
            IIC = IROW(LNC+1)
            do I=1,NVIRC
              RSUM = DDOT_(NVIRA,A(INN),1,A(INN),1)
              RSUM = COPI*RSUM
              TSUM = TSUM+RSUM
              IIC = IIC+LNC+I
              FC(IIC) = FC(IIC)+RSUM
              INN = INN+NVIRA
            end do
          end if
        end do
        if (IDENS == 0) cycle
        TSUM = TSUM*Half
      else
        ! DOUBLET-DOUBLET INTERACTIONS
        if (NVIR(MYL) == 0) cycle
        if (IDENS /= 1) then
          call FMMM(F(IPOF(MYL)+1),C(INMY),A,NVIR(MYL),1,NVIR(MYL))
          S(INMY:INMY+NVIR(MYL)-1) = S(INMY:INMY+NVIR(MYL)-1)+FACS*A(1:NVIR(MYL))
          W(INMY:INMY+NVIR(MYL)-1) = W(INMY:INMY+NVIR(MYL)-1)+FACW*A(1:NVIR(MYL))
          cycle
        else
          call FMUL2(C(INMY),C(INMY),A,NVIR(MYL),NVIR(MYL),1)
          IPF = IPOF(MYL)+1
          IIN = IPOF(MYL+1)-IPOF(MYL)
          ENPQ = (One-THET(INDA,INDA)*Half)*(ENP(INDA)+ENP(INDA)-One)+THET(INDA,INDA)*Half
          COPI = One/ENPQ
          F(IPF:IPF+IIN-1) = F(IPF:IPF+IIN-1)+COPI*A(1:IIN)
          NVIRA = NVIR(MYL)
          LNA = LN+NSYS(MYL)
          IIA = IROW(LNA+1)
          TSUM = Zero
          do I=1,NVIRA
            RSUM = COPI*C(INMY)*C(INMY)
            INMY = INMY+1
            TSUM = TSUM+RSUM
            IIA = IIA+LNA+I
            FC(IIA) = FC(IIA)+RSUM
          end do
        end if
      end if
    end if
    IJ = 0
    do I=1,LN
      IJ = IJ+I
      FC(IJ) = FC(IJ)+IOC(I)*TSUM
    end do
  end do
  ITURN = 1
  if (IDENS /= 1) exit
end do
call MDSQ2(C,S,W,MUL,INDX,JSY,NDIAG,INUM,IRC(3),LSYM,NVIRT,SQ2)
if (IPRINT >= 15) then
  write(u6,'(A,/,(10F12.6))') ' S,AB',(S(I),I=1,JSC(4))
  write(u6,'(A,/,(10F12.6))') ' W,AB',(W(I),I=1,JSC(4))
  if (IDENS == 1) write(u6,'(A,/,(10F12.6))') ' FC,AB',(FC(I),I=1,NOB2)
end if

return

310 format(/,6X,'TRACE OF DENSITY MATRIX',F16.8)

end subroutine MAB
