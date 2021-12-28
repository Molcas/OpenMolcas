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

subroutine INT4(I,J,L,IT1,IT2,II,IID,JJ,JJD,JTYP,ITAI,L0,L1,L2,L3)
! I < J == K < L

use guga_global, only: COUP, IJ, IWAY, IX, J1, J2
use Definitions, only: wp, iwp

#include "intent.fh"

implicit none
integer(kind=iwp), intent(in) :: I, J, L, IT1, IT2, II, IID, JJ, JJD, JTYP, L0(*), L1(*), L2(*), L3(*)
integer(kind=iwp), intent(_OUT_) :: ITAI(*)
integer(kind=iwp) :: ISTOP, ITAIL, ITURN, ITYP, KM, LJ, LJM, LJS
logical(kind=iwp) :: first1, first2, skip

ITYP = 0
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
    end if
    KM = KM-1
    IWAY(KM) = 1
    first1 = KM == J
    loop_2: do
      if (first1) then
        first1 = .false.
        ITURN = 0
        IWAY(J) = 1
        loop_3: do
          KM = J
          if (ITURN == 0) then
            call LOOP12(KM,ISTOP,IT1,IT2)
          else if (ITURN == 1) then
            call LOOP9(KM,ISTOP,IT1,IT2)
          end if
          if (ISTOP /= 0) then
            if (ITURN /= 1) then
              ITURN = 1
              IWAY(J) = 1
              cycle loop_3
            else
              KM = KM+1
              if (KM == L) cycle loop_1
              cycle loop_2
            end if
          end if
          loop_4: do
            KM = KM-1
            first2 = .false.
            if (KM == 0) then
              call COMP1(LJ,ITYP,L,IT2,II,IID,JJ,JJD,JTYP,ITAI)
              KM = 1
              if (KM == J) cycle loop_3
            else
              IWAY(KM) = 1
              if (KM == I) first2 = .true.
            end if
            loop_5: do
              if (first2) then
                IWAY(I) = 1
                do
                  KM = I
                  if (ITURN == 0) then
                    call LOOP3(KM,ISTOP,IT1,IT2)
                  else if (ITURN == 1) then
                    call LOOP4(KM,ISTOP,IT1,IT2)
                  end if
                  if (ISTOP == 1) exit
                  if (abs(COUP(I)) < 1.0e-6_wp) cycle
                  call COMP(I,LJ,ITYP,I,IT1,IT2)
                end do
                first2 = .false.
              else
                if (ITURN == 0) then
                  call LOOP5(KM,ISTOP,IT1,IT2)
                else if (ITURN == 1) then
                  call LOOP6(KM,ISTOP,IT1,IT2)
                end if
                if (ISTOP == 0) cycle loop_4
              end if
              KM = KM+1
              if (KM == J) cycle loop_3
            end do loop_5
            if (.true.) exit loop_4
          end do loop_4
          if (.true.) exit loop_3
        end do loop_3
        exit loop_2
      else
        call LOOP5(KM,ISTOP,IT1,IT2)
        if (ISTOP == 0) then
          skip = .true.
          cycle loop_1
        end if
        KM = KM+1
        if (KM == L) cycle loop_1
      end if
    end do loop_2
    if (.true.) exit loop_1
  end do loop_1
end do

return

end subroutine INT4
