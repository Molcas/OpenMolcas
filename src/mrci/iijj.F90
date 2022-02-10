!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!***********************************************************************

subroutine IIJJ(ICSPCK,INTSYM,HDIAG,FC,FIIJJ,FIJIJ)

implicit real*8(A-H,O-Z)
#include "SysDef.fh"
#include "mrci.fh"
dimension ICSPCK(*), INTSYM(*), HDIAG(*), FC(*), FIIJJ(*), FIJIJ(*)
dimension IOC(55)
!Statement functions
JO(L) = ICUNP(ICSPCK,L)
JSYM(L) = JSUNP(INTSYM,L)

IAD27 = 0
II1 = 0
ILIM = 4
if (IFIRST /= 0) ILIM = 2
IRL = IRC(ILIM)
do IR=1,IRL
  do I=1,LN
    II1 = II1+1
    JOJ = JO(II1)
    if (JOJ > 1) JOJ = JOJ-1
    IOC(I) = JOJ
  end do
  NSS = MUL(JSYM(IR),LSYM)
  SUM = 0.0d00
  do I=1,LN
    IJ = IROW(I)
    if (IOC(I) == 0) GO TO 111
    do J=1,I-1
      IJ = IJ+1
      if (IOC(J) == 0) GO TO 113
      SUM = SUM+IOC(I)*(IOC(J)*FIIJJ(IJ)-FIJIJ(IJ))
113   continue
    end do
    IJ = IJ+1
    SUM = SUM+(IOC(I)-1)*FIIJJ(IJ)+IOC(I)*FC(IJ)
111 continue
  end do
  if (IR > IRC(1)) GO TO 120
  HDIAG(IR) = SUM
  if (IR /= IRC(1)) GO TO 100
  call dDAFILE(Lu_27,1,HDIAG,IRC(1),IAD27)
  GO TO 100
120 IND = 0
  if (IR > IRC(2)) goto 130
  NA1 = NVIRP(NSS)+1
  NA2 = NVIRP(NSS)+NVIR(NSS)
  if (NA2 < NA1) GO TO 100
  do NA=NA1,NA2
    IND = IND+1
    IA = IROW(LN+NA)
    SUM1 = SUM+FC(IA+LN+NA)
    do I=1,LN
      if (IOC(I) == 0) GO TO 122
      SUM1 = SUM1+IOC(I)*FIIJJ(IA+I)-FIJIJ(IA+I)
122   continue
    end do
    HDIAG(IND) = SUM1
  end do
  call dDAFILE(Lu_27,1,HDIAG,IND,IAD27)
  GO TO 100
130 IND = 0
  do NA=1,NVIRT
    NSA = MUL(NSS,NSM(LN+NA))
    NB1 = NVIRP(NSA)+1
    NB2 = NVIRP(NSA)+NVIR(NSA)
    if (NB2 > NA) NB2 = NA
    if (NB2 < NB1) GO TO 141
    IA = IROW(LN+NA)
    IAV = IA+LN
    do NB=NB1,NB2
      IND = IND+1
      IB = IROW(LN+NB)
      IBV = IB+LN
      TERM = SUM+FIIJJ(IAV+NB)+FC(IAV+NA)+FC(IBV+NB)
      if (IR <= IRC(3)) then
        SUM1 = TERM-FIJIJ(IAV+NB)
      else
        SUM1 = TERM+FIJIJ(IAV+NB)
      end if
      do I=1,LN
        if (IOC(I) == 0) GO TO 143
        TERM = IOC(I)*(FIIJJ(IA+I)+FIIJJ(IB+I))-FIJIJ(IA+I)-FIJIJ(IB+I)
        SUM1 = SUM1+TERM
143     continue
      end do
      HDIAG(IND) = SUM1
    end do
141 continue
  end do
  if (IND > 0) call dDAFILE(Lu_27,1,HDIAG,IND,IAD27)
100 continue
end do

return

end subroutine IIJJ
