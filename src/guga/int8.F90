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

subroutine INT8(I,J,L,IT1,IT2,II,IID,JJT,JJD,JTYP,ITAI,L0,L1,L2,L3)
! K == L , I < J

use guga_global, only: COUP, IJ, IWAY, IX, J1, J2
use Constants, only: One, Two
use Definitions, only: wp, iwp

#include "intent.fh"

implicit none
integer(kind=iwp), intent(in) :: I, J, L, IT1, IT2, II, IID, JJT, JJD, JTYP, L0(*), L1(*), L2(*), L3(*)
integer(kind=iwp), intent(_OUT_) :: ITAI(*)
integer(kind=iwp) :: ISTOP, ITAIL, ITYP, JJ, JJM, JJS, KM
real(kind=wp) :: FAC
logical(kind=iwp) :: first, skip

ITYP = 0
if ((L < I) .or. (L > J)) ITYP = 1
if (I == 0) ITYP = 0
if ((I == 0) .and. (L > J)) ITYP = 1
FAC = One
JJS = IJ(J+1)+1
JJM = IJ(J)
do JJ=JJS,JJM
  ITAIL = IX(IT2+JJ)
  if (IT1 /= IT2) call TAIL(J,JJ,ITAI,ITAIL,L0,L1,L2,L3,IT1,IT2)
  IWAY(J) = 1
  skip = .false.
  loop_1: do
    if (skip) then
      skip = .false.
    else
      KM = J
      J2(KM+1) = JJ
      J1(KM+1) = JJ
      call LOOP1(KM,ISTOP,IT1,IT2)
      if (ISTOP == 1) exit loop_1
    end if
    KM = KM-1
    first = .false.
    if (KM == 0) then
      COUP(1) = FAC*COUP(1)
      call COMP1(JJ,ITYP,L,IT2,II,IID,JJT,JJD,JTYP,ITAI)
      KM = 1
      if (KM == J) cycle loop_1
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
          if (abs(COUP(I)) < 1.0e-6_wp) cycle
          COUP(I) = FAC*COUP(I)
          call COMP(I,JJ,ITYP,L,IT1,IT2)
        end do
        first = .false.
      else
        call LOOP5(KM,ISTOP,IT1,IT2)
        if (ISTOP /= 1) then
          if (KM == L) then
            if (IWAY(KM) == 2) cycle
            FAC = One
            if (IWAY(KM) == 5) FAC = Two
          end if
          skip = .true.
          cycle loop_1
        end if
      end if
      KM = KM+1
      if (KM == J) cycle loop_1
    end do
    if (.true.) exit loop_1
  end do loop_1
end do

return

end subroutine INT8
