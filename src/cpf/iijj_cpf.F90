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

subroutine IIJJ_CPF(ICASE,JSY,HDIAG,FC,FIJ,FJI)

implicit real*8(A-H,O-Z)
#include "SysDef.fh"
#include "cpfmcpf.fh"
#include "files_cpf.fh"
dimension JSY(*), HDIAG(*), FC(*), FIJ(*), FJI(*)
dimension ICASE(*)
dimension IOC(55)
! Statement functions
JO(L) = ICUNP(ICASE,L)
JSYM(L) = JSUNP_CPF(JSY,L)

IAD27 = 0
ILIM = 4
if (IFIRST /= 0) ILIM = 2
IRL = IRC(ILIM)

do IR=1,IRL
  do I=1,LN
    JOJ = JO(I+LN*(IR-1))
    IOC(I) = (JOJ+1)/2
  end do
  NSS = MUL(JSYM(IR),LSYM)

  SUM = 0.0d0
  do I=1,LN
    if (IOC(I) /= 0) then
      do J=1,I-1
        IJ = (I*(I-1))/2+J
        if (IOC(J) /= 0) then
          TERM = IOC(I)*(IOC(J)*FIJ(IJ)-FJI(IJ))
          SUM = SUM+TERM
        end if
      end do
      II = (I*(I+1))/2
      TERM = (IOC(I)-1)*FIJ(II)+IOC(I)*FC(II)
      SUM = SUM+TERM
    end if
  end do

  if (IR > IRC(1)) GO TO 120
  ! IR=1..IRC(1), HDIAG(IR)=SUM
  HDIAG(IR) = SUM
  if (IR == IRC(1)) call dDAFILE(Lu_27,1,HDIAG,IRC(1),IAD27)
  GO TO 100

120 continue
  if (IR > IRC(2)) GO TO 130
  ! IR=IRC(1)+1 ... IRC(2)
  IND = 0
  NA1 = NSYS(NSS)+1
  NA2 = NSYS(NSS+1)
  if (NA2 < NA1) GO TO 100
  do NA=NA1,NA2
    IND = IND+1
    IA = IROW(LN+NA)
    SUM1 = SUM+FC(IA+LN+NA)
    do I=1,LN
      if (IOC(I) /= 0) SUM1 = SUM1+IOC(I)*FIJ(IA+I)-FJI(IA+I)
    end do
    HDIAG(IND) = SUM1
  end do
  call dDAFILE(Lu_27,1,HDIAG,IND,IAD27)
  GO TO 100

130 continue
  ! IR=IRC(2)+1 ... IRC(ILIM)
  IND = 0
  do NA=1,NVIRT
    NSA = MUL(NSS,NSM(LN+NA))
    NB1 = NSYS(NSA)+1
    NB2 = NSYS(NSA+1)
    if (NB2 > NA) NB2 = NA
    if (NB2 < NB1) GO TO 141
    IA = IROW(LN+NA)
    IAV = IA+LN
    do NB=NB1,NB2
      IND = IND+1
      IB = IROW(LN+NB)
      IBV = IB+LN
      TERM = SUM+FIJ(IAV+NB)+FC(IAV+NA)+FC(IBV+NB)
      if (IR <= IRC(3)) SUM1 = TERM-FJI(IAV+NB)
      if (IR > IRC(3)) SUM1 = TERM+FJI(IAV+NB)
      do I=1,LN
        if (IOC(I) /= 0) then
          TERM = IOC(I)*(FIJ(IA+I)+FIJ(IB+I))-FJI(IA+I)-FJI(IB+I)
          SUM1 = SUM1+TERM
        end if
      end do
      HDIAG(IND) = SUM1
    end do
141 continue
  end do
  if (IND > 0) call dDAFILE(Lu_27,1,HDIAG,IND,IAD27)
100 continue
end do

return

end subroutine IIJJ_CPF
