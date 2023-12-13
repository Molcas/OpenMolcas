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

subroutine INT9(I,J,L,IT1,IT2,II,IID,JJ,JJD,JTYP,ITAI,L0,L1,L2,L3)
! I < L . CASE 1 J=K=L , CASE 2 I=J=K .

use guga_global, only: COUP, IJ, IWAY, IX, J1, J2
use Constants, only: Zero, One
use Definitions, only: wp, iwp

#include "intent.fh"

implicit none
integer(kind=iwp), intent(in) :: I, J, L, IT1, IT2, II, IID, JJ, JJD, JTYP, L0(*), L1(*), L2(*), L3(*)
integer(kind=iwp), intent(_OUT_) :: ITAI(*)
integer(kind=iwp) :: ISTOP, ITAIL, ITYP, KM, LJ, LJM, LJS
real(kind=wp) :: FAC
logical(kind=iwp) :: first, skip

ITYP = 0
FAC = One
LJS = IJ(L+1)+1
LJM = IJ(L)
do LJ=LJS,LJM
  ITAIL = IX(IT2+LJ)
  if (IT1 /= IT2) call TAIL(L,LJ,ITAI,ITAIL,L0,L1,L2,L3,IT1,IT2)
  IWAY(L) = 1
  skip = .false.
  loop_1: do
    if (skip) then
      skip = .false.
    else
      KM = L
      J2(KM+1) = LJ
      J1(KM+1) = LJ
      call LOOP1(KM,ISTOP,IT1,IT2)
      if (ISTOP == 1) exit loop_1
      if (I /= J) then
        FAC = One
        if ((IWAY(KM) == 2) .or. (IWAY(KM) == 4)) FAC = Zero
      end if
    end if
    KM = KM-1
    first = .false.
    if (KM == 0) then
      if (FAC /= Zero) call COMP1(LJ,ITYP,L,IT2,II,IID,JJ,JJD,JTYP,ITAI)
      KM = 1
      if (KM == L) cycle loop_1
    else
      IWAY(KM) = 1
      if (KM == I) first = .true.
    end if
    do
      if (first) then
        IWAY(I) = 1
        do
          KM = I
          call LOOP3(KM,ISTOP,IT1,IT2)
          if (ISTOP == 1) exit
          if (I == J) then
            FAC = One
            if (IWAY(KM) == 2) FAC = Zero
          end if
          if (FAC == Zero) cycle
          if (abs(COUP(I)) < 1.0e-6_wp) cycle
          call COMP(I,LJ,ITYP,I,IT1,IT2)
        end do
        first = .false.
      else
        call LOOP5(KM,ISTOP,IT1,IT2)
        if (ISTOP == 0) then
          skip = .true.
          cycle loop_1
        end if
      end if
      KM = KM+1
      if (KM == L) cycle loop_1
    end do
    if (.true.) exit loop_1
  end do loop_1
end do

return

end subroutine INT9
