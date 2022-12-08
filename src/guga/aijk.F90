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
!               2021, Ignacio Fdez. Galvan                             *
!***********************************************************************
! 2021: Remove GOTOs

subroutine AIJK(ITAI,L0,L1,L2,L3)

use guga_global, only: IADD10, ICH, IOUT, IRC, JRC, LN, Lu_10, MXVERT, NBUF, NMAT
use guga_util_global, only: COP, ICOP1, nCOP
use Definitions, only: iwp, u6

#include "intent.fh"

implicit none
integer(kind=iwp), intent(_OUT_) :: ITAI(*)
integer(kind=iwp), intent(in) :: L0(*), L1(*), L2(*), L3(*)
integer(kind=iwp) :: I, II, IID, IND, IT1, IT2, ITT1, ITT2, ITURN, J, JJ, JJD, K, L, NI, NJ, NK

IOUT = 0
NMAT = 0
L = 0
do NI=1,LN
  do NJ=1,NI
    I = ICH(NI)
    J = ICH(NJ)
    if (I <= J) then
      I = ICH(NJ)
      J = ICH(NI)
    end if
    do NK=1,LN
      IOUT = IOUT+1
      ICOP1(IOUT) = 0
      if (IOUT >= NBUF) then
        ICOP1(nCOP+1) = NBUF
        call dDAFILE(Lu_10,1,COP,NCOP,IADD10)
        call iDAFILE(Lu_10,1,iCOP1,NCOP+1,IADD10)
        NMAT = NMAT+NBUF
        IOUT = 0
      end if
      K = ICH(NK)
      IOUT = IOUT+1
      !IND = I+2**10*J
      IND = ior(I,ishft(J,10))
      !ICOP1(IOUT) = IND+2**20*K
      ICOP1(IOUT) = ior(IND,ishft(K,20))
      if (IOUT >= NBUF) then
        ICOP1(nCOP+1) = NBUF
        call dDAFILE(Lu_10,1,COP,NCOP,IADD10)
        call iDAFILE(Lu_10,1,iCOP1,NCOP+1,IADD10)
        NMAT = NMAT+NBUF
        IOUT = 0
      end if
      ! DOUBLET-VALENCE INTERACTIONS
      ITURN = 1
      ITT1 = 1
      ITT2 = 0
      do
        IT1 = ITT1*MXVERT
        IT2 = ITT2*MXVERT
        if (ITT1 == 0) then
          JJ = 0
          JJD = 0
        else
          JJ = IRC(ITT1)
          JJD = JRC(ITT1)
        end if
        if (ITT2 == 0) then
          II = 0
          IID = 0
        else
          II = IRC(ITT2)
          IID = JRC(ITT2)
        end if
        if (I == K) then
          if (J /= K) then
            call INT61(L,J,I,IT1,IT2,II,IID,JJ,JJD,ITURN,ITAI,L0,L1,L2,L3)
            call INT62(L,J,I,IT1,IT2,II,IID,JJ,JJD,ITURN,ITAI,L0,L1,L2,L3)
          else
            call INT9(L,K,I,IT1,IT2,II,IID,JJ,JJD,ITURN,ITAI,L0,L1,L2,L3)
          end if
        else if (J == K) then
          call INT4(L,K,I,IT1,IT2,II,IID,JJ,JJD,ITURN,ITAI,L0,L1,L2,L3)
        else if (I == J) then
          call INT8(L,K,I,IT1,IT2,II,IID,JJ,JJD,ITURN,ITAI,L0,L1,L2,L3)
        else if (K >= I) then
          call INT3(J,I,L,K,IT1,IT2,II,IID,JJ,JJD,ITURN,ITAI,L0,L1,L2,L3)
        else if (K >= J) then
          call INT2(L,K,J,I,IT1,IT2,II,IID,JJ,JJD,ITURN,ITAI,L0,L1,L2,L3)
        else
          call INT1(L,K,J,I,IT1,IT2,II,IID,JJ,JJD,ITURN,ITAI,L0,L1,L2,L3)
        end if
        select case (ITURN)
          case default !(1)
            ! TRIPLET-DOUBLET INTERACTIONS
            ITURN = 2
            ITT1 = 2
            ITT2 = 1
          case (2)
            ! SINGLET-DOUBLET INTERACTIONS
            ITURN = 3
            ITT1 = 3
            ITT2 = 1
          case (3)
            ! VALENCE-DOUBLET INTERACTIONS
            ITURN = 4
            ITT1 = 0
            ITT2 = 1
          case (4)
            ! DOUBLET-TRIPLET INTERACTIONS
            ITURN = 5
            ITT1 = 1
            ITT2 = 2
          case (5)
            ! DOUBLET-SINGLET INTERACTIONS
            ITURN = 6
            ITT1 = 1
            ITT2 = 3
          case (6)
            exit
        end select
      end do
    end do
  end do
end do
ICOP1(nCOP+1) = IOUT
call dDAFILE(Lu_10,1,COP,NCOP,IADD10)
call iDAFILE(Lu_10,1,iCOP1,NCOP+1,IADD10)
NMAT = NMAT+IOUT
ICOP1(nCOP+1) = -1
call dDAFILE(Lu_10,1,COP,NCOP,IADD10)
call iDAFILE(Lu_10,1,iCOP1,NCOP+1,IADD10)
write(u6,600) NMAT

return

600 format(/6X,'COEFFICIENTS FOR AIJK',I9)

end subroutine AIJK
