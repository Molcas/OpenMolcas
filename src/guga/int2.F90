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

subroutine INT2(I,J,K,L,IT1,IT2,II,IID,JJ,JJD,JTYP,ITAI,L0,L1,L2,L3)
! I < K < J < L

use guga_global, only: COUP, IJ, IVF0, IWAY, IX, J1, J2, JM, JM1
use Definitions, only: wp, iwp

#include "intent.fh"

implicit none
integer(kind=iwp), intent(in) :: I, J, K, L, IT1, IT2, II, IID, JJ, JJD, JTYP, L0(*), L1(*), L2(*), L3(*)
integer(kind=iwp), intent(_OUT_) :: ITAI(*)
integer(kind=iwp) :: ISTOP, ITAIL, ITURN, ITYP, KM, LJ, LJM, LJS
logical(kind=iwp) :: first1, first2, first3

ITYP = 0
LJS = IJ(L+1)+1
LJM = IJ(L)
do LJ=LJS,LJM
  ITAIL = IX(IT2+LJ)
  if (IT1 /= IT2) call TAIL(L,LJ,ITAI,ITAIL,L0,L1,L2,L3,IT1,IT2)
  IWAY(L) = 1
  loop_1: do
    KM = L
    J2(KM+1) = LJ
    J1(KM+1) = LJ
    call LOOP1(KM,ISTOP,IT1,IT2)
    if (ISTOP == 1) exit loop_1
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
          JM(KM) = IVF0+1
          JM1(KM) = IVF0+1
          if (ITURN == 0) then
            call LOOP10(KM,ISTOP,IT1,IT2)
          else if (ITURN == 1) then
            call LOOP11(KM,ISTOP,IT1,IT2)
          end if
          first2 = .false.
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
          else
            KM = KM-1
            IWAY(KM) = 1
            if (KM == K) first2 = .true.
          end if
          loop_4: do
            if (first2) then
              first2 = .false.
              IWAY(K) = 1
              loop_5: do
                KM = K
                if (ITURN == 0) then
                  call LOOP16(KM,ISTOP,IT1,IT2)
                else if (ITURN == 1) then
                  call LOOP20(KM,ISTOP,IT1,IT2)
                end if
                if (ISTOP /= 0) then
                  KM = KM+1
                  if (KM == J) cycle loop_3
                  cycle loop_4
                end if
                loop_6: do
                  KM = KM-1
                  first3 = .false.
                  if (KM == 0) then
                    call COMP1(LJ,ITYP,L,IT2,II,IID,JJ,JJD,JTYP,ITAI)
                    KM = 1
                    if (KM == K) cycle loop_5
                  else
                    IWAY(KM) = 1
                    if (KM == I) first3 = .true.
                  end if
                  do
                    if (first3) then
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
                      first3 = .false.
                    else
                      if (ITURN == 0) then
                        call LOOP5(KM,ISTOP,IT1,IT2)
                      else if (ITURN == 1) then
                        call LOOP6(KM,ISTOP,IT1,IT2)
                      end if
                      if (ISTOP == 0) cycle loop_6
                    end if
                    KM = KM+1
                    if (KM == K) cycle loop_5
                  end do
                  if (.true.) exit loop_6
                end do loop_6
                if (.true.) exit loop_5
              end do loop_5
              exit loop_4
            else
              JM(KM) = IVF0+1
              JM1(KM) = IVF0+1
              if (ITURN == 0) then
                call LOOP17(KM,ISTOP,IT1,IT2)
              else if (ITURN == 1) then
                call LOOP21(KM,ISTOP,IT1,IT2)
              end if
              if (ISTOP == 0) then
                KM = KM-1
                IWAY(KM) = 1
                if (KM == K) first2 = .true.
              else
                KM = KM+1
                if (KM == J) cycle loop_3
              end if
            end if
          end do loop_4
          if (.true.) exit loop_3
        end do loop_3
        exit loop_2
      else
        call LOOP5(KM,ISTOP,IT1,IT2)
        if (ISTOP == 0) then
          KM = KM-1
          IWAY(KM) = 1
          if (KM == J) first1 = .true.
        else
          KM = KM+1
          if (KM == L) cycle loop_1
        end if
      end if
    end do loop_2
    if (.true.) exit loop_1
  end do loop_1
end do

return

end subroutine INT2
