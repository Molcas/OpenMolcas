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

use cpf_global, only: ILIM, IRC, IROW, LN, LSYM, Lu_27, NSM, NSYS, NVIRT
use Symmetry_Info, only: Mul
use Constants, only: Zero
use Definitions, only: wp, iwp

#include "intent.fh"

implicit none
integer(kind=iwp), intent(in) :: ICASE(*), JSY(*)
real(kind=wp), intent(_OUT_) :: HDIAG(*)
real(kind=wp), intent(in) :: FC(*), FIJ(*), FJI(*)
integer(kind=iwp) :: I, IA, IAD27, IAV, IB, IBV, II, IJ, IND, IOC(55), IR, IRL, J, JOJ, NA, NA1, NA2, NB, NB1, NB2, NSA, NSS
real(kind=wp) :: SUM1, SUM2, TERM
integer(kind=iwp), external :: ICUNP, JSUNP

IAD27 = 0
IRL = IRC(ILIM)

do IR=1,IRL
  do I=1,LN
    JOJ = ICUNP(ICASE,I+LN*(IR-1))
    IOC(I) = (JOJ+1)/2
  end do
  NSS = MUL(JSUNP(JSY,IR),LSYM)

  SUM1 = Zero
  do I=1,LN
    if (IOC(I) /= 0) then
      do J=1,I-1
        IJ = (I*(I-1))/2+J
        if (IOC(J) /= 0) then
          TERM = IOC(I)*(IOC(J)*FIJ(IJ)-FJI(IJ))
          SUM1 = SUM1+TERM
        end if
      end do
      II = (I*(I+1))/2
      TERM = (IOC(I)-1)*FIJ(II)+IOC(I)*FC(II)
      SUM1 = SUM1+TERM
    end if
  end do

  if (IR <= IRC(1)) then

    ! IR=1..IRC(1), HDIAG(IR)=SUM1
    HDIAG(IR) = SUM1
    if (IR == IRC(1)) call dDAFILE(Lu_27,1,HDIAG,IRC(1),IAD27)

  else if (IR <= IRC(2)) then

    ! IR=IRC(1)+1 ... IRC(2)
    IND = 0
    NA1 = NSYS(NSS)+1
    NA2 = NSYS(NSS+1)
    if (NA2 >= NA1) then
      do NA=NA1,NA2
        IND = IND+1
        IA = IROW(LN+NA)
        SUM2 = SUM1+FC(IA+LN+NA)
        do I=1,LN
          if (IOC(I) /= 0) SUM2 = SUM2+IOC(I)*FIJ(IA+I)-FJI(IA+I)
        end do
        HDIAG(IND) = SUM2
      end do
      call dDAFILE(Lu_27,1,HDIAG,IND,IAD27)
    end if

  else

    ! IR=IRC(2)+1 ... IRC(ILIM)
    IND = 0
    do NA=1,NVIRT
      NSA = MUL(NSS,NSM(LN+NA))
      NB1 = NSYS(NSA)+1
      NB2 = NSYS(NSA+1)
      if (NB2 > NA) NB2 = NA
      if (NB2 >= NB1) then
        IA = IROW(LN+NA)
        IAV = IA+LN
        do NB=NB1,NB2
          IND = IND+1
          IB = IROW(LN+NB)
          IBV = IB+LN
          TERM = SUM1+FIJ(IAV+NB)+FC(IAV+NA)+FC(IBV+NB)
          if (IR <= IRC(3)) SUM2 = TERM-FJI(IAV+NB)
          if (IR > IRC(3)) SUM2 = TERM+FJI(IAV+NB)
          do I=1,LN
            if (IOC(I) /= 0) then
              TERM = IOC(I)*(FIJ(IA+I)+FIJ(IB+I))-FJI(IA+I)-FJI(IB+I)
              SUM2 = SUM2+TERM
            end if
          end do
          HDIAG(IND) = SUM2
        end do
      end if
    end do
    if (IND > 0) call dDAFILE(Lu_27,1,HDIAG,IND,IAD27)

  end if
end do

return

end subroutine IIJJ_CPF
