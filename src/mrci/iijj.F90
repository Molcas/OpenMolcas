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

use mrci_global, only: IFIRST, IRC, IROW, LN, LSYM, Lu_27, NSM, NVIR, NVIRP, NVIRT
use Symmetry_Info, only: Mul
use Constants, only: Zero
use Definitions, only: wp, iwp

#include "intent.fh"

implicit none
integer(kind=iwp), intent(in) :: ICSPCK(*), INTSYM(*)
real(kind=wp), intent(_OUT_) :: HDIAG(*)
real(kind=wp), intent(in) :: FC(*), FIIJJ(*), FIJIJ(*)
integer(kind=iwp) :: I, IA, IAD27, IAV, IB, IBV, II1, IJ, ILIM, IND, IOC(55), IR, IRL, J, JOJ, NA, NA1, NA2, NB, NB1, NB2, NSA, NSS
real(kind=wp) :: SUM1, SUM2, TERM
integer(kind=iwp), external :: ICUNP, JSUNP

IAD27 = 0
II1 = 0
ILIM = 4
if (IFIRST /= 0) ILIM = 2
IRL = IRC(ILIM)
do IR=1,IRL
  do I=1,LN
    II1 = II1+1
    JOJ = ICUNP(ICSPCK,II1)
    if (JOJ > 1) JOJ = JOJ-1
    IOC(I) = JOJ
  end do
  NSS = MUL(JSUNP(INTSYM,IR),LSYM)
  SUM1 = Zero
  do I=1,LN
    IJ = IROW(I)
    if (IOC(I) == 0) cycle
    do J=1,I-1
      IJ = IJ+1
      if (IOC(J) /= 0) SUM1 = SUM1+IOC(I)*(IOC(J)*FIIJJ(IJ)-FIJIJ(IJ))
    end do
    IJ = IJ+1
    SUM1 = SUM1+(IOC(I)-1)*FIIJJ(IJ)+IOC(I)*FC(IJ)
  end do
  if (IR <= IRC(1)) then
    HDIAG(IR) = SUM1
    if (IR == IRC(1)) call dDAFILE(Lu_27,1,HDIAG,IRC(1),IAD27)
  else
    IND = 0
    if (IR <= IRC(2)) then
      NA1 = NVIRP(NSS)+1
      NA2 = NVIRP(NSS)+NVIR(NSS)
      do NA=NA1,NA2
        IND = IND+1
        IA = IROW(LN+NA)
        SUM2 = SUM1+FC(IA+LN+NA)
        do I=1,LN
          if (IOC(I) /= 0) SUM2 = SUM2+IOC(I)*FIIJJ(IA+I)-FIJIJ(IA+I)
        end do
        HDIAG(IND) = SUM2
      end do
      call dDAFILE(Lu_27,1,HDIAG,IND,IAD27)
    else
      IND = 0
      do NA=1,NVIRT
        NSA = MUL(NSS,NSM(LN+NA))
        NB1 = NVIRP(NSA)+1
        NB2 = NVIRP(NSA)+NVIR(NSA)
        if (NB2 > NA) NB2 = NA
        if (NB2 < NB1) cycle
        IA = IROW(LN+NA)
        IAV = IA+LN
        do NB=NB1,NB2
          IND = IND+1
          IB = IROW(LN+NB)
          IBV = IB+LN
          TERM = SUM1+FIIJJ(IAV+NB)+FC(IAV+NA)+FC(IBV+NB)
          if (IR <= IRC(3)) then
            SUM2 = TERM-FIJIJ(IAV+NB)
          else
            SUM2 = TERM+FIJIJ(IAV+NB)
          end if
          do I=1,LN
            if (IOC(I) /= 0) then
              TERM = IOC(I)*(FIIJJ(IA+I)+FIIJJ(IB+I))-FIJIJ(IA+I)-FIJIJ(IB+I)
              SUM2 = SUM2+TERM
            end if
          end do
          HDIAG(IND) = SUM2
        end do
      end do
      if (IND > 0) call dDAFILE(Lu_27,1,HDIAG,IND,IAD27)
    end if
  end if
end do

return

end subroutine IIJJ
