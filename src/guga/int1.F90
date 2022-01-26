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

subroutine INT1(I,J,K,L,IT1,IT2,II,IID,JJ,JJD,JTYP,ITAI,L0,L1,L2,L3)
! I < J < K < L

use guga_global, only: COUP, ICOUP, ICOUP1, IJ, IWAY, IX, J1, J2
use Definitions, only: wp, iwp

#include "intent.fh"

implicit none
integer(kind=iwp), intent(in) :: I, J, K, L, IT1, IT2, II, IID, JJ, JJD, JTYP, L0(*), L1(*), L2(*), L3(*)
integer(kind=iwp), intent(_OUT_) :: ITAI(*)
integer(kind=iwp) :: ICP1I, ICPI, ISTOP, ITAIL, ITYP, KM, LJ, LJM, LJS
logical(kind=iwp) :: first1, first2, first3, skip

LJS = IJ(L+1)+1
LJM = IJ(L)
ITYP = 0
do LJ=LJS,LJM
  ITAIL = IX(IT2+LJ)
  if (IT1 /= IT2) call TAIL(L,LJ,ITAI,ITAIL,L0,L1,L2,L3,IT1,IT2)
  IWAY(L) = 1
  loop_1: do
    KM = L
    J2(KM+1) = LJ
    J1(KM+1) = LJ
    call LOOP1(L,ISTOP,IT1,IT2)
    if (ISTOP == 1) exit loop_1
    first1 = .true.
    loop_2: do
      if (first1) then
        KM = KM-1
        IWAY(KM) = 1
        first1 = .false.
        if (KM == K) then
          IWAY(K) = 1
          skip = .false.
          loop_3: do
            if (.not. skip) then
              KM = K
              call LOOP3(KM,ISTOP,IT1,IT2)
              if (ISTOP /= 0) then
                KM = KM+1
                if (KM == L) cycle loop_1
                cycle loop_2
              end if
            end if
            skip = .false.
            KM = KM-1
            IWAY(KM) = 1
            first2 = KM == J
            loop_4: do
              if (first2) then
                KM = J
                if (IT1 >= IT2) then
                  call LOOP1(KM,ISTOP,IT1,IT2)
                else
                  call LOOP2(KM,ISTOP,IT1,IT2)
                end if
                if (ISTOP /= 0) then
                  KM = KM+1
                  if (KM == K) cycle loop_3
                  first2 = .false.
                  cycle loop_4
                end if
                loop_5: do
                  KM = KM-1
                  first3 = .false.
                  if (KM == 0) then
                    COUP(1) = COUP(1)*COUP(K)
                    if (abs(COUP(1)) >= 1.0e-6_wp) then
                      ICOUP(1) = ICOUP(J+1)+ICOUP(1)
                      ICOUP1(1) = ICOUP1(J+1)+ICOUP1(1)
                      call COMP1(LJ,ITYP,L,IT2,II,IID,JJ,JJD,JTYP,ITAI)
                    end if
                    KM = 1
                    if (KM == J) cycle loop_4
                  else
                    IWAY(KM) = 1
                    if (KM == I) first3 = .true.
                  end if
                  do
                    if (first3) then
                      IWAY(I) = 1
                      do
                        KM = I
                        if (IT1 >= IT2) then
                          call LOOP3(KM,ISTOP,IT1,IT2)
                        else
                          call LOOP4(KM,ISTOP,IT1,IT2)
                        end if
                        if (ISTOP == 1) exit
                        COUP(I) = COUP(I)*COUP(K)
                        if (abs(COUP(I)) < 1.0e-6_wp) cycle
                        ICPI = ICOUP(I)
                        ICP1I = ICOUP1(I)
                        ICOUP(I) = ICOUP(J+1)+ICPI
                        ICOUP1(I) = ICOUP1(J+1)+ICP1I
                        call COMP(I,LJ,ITYP,I,IT1,IT2)
                        ICOUP(I) = ICOUP(J+1)+ICP1I
                        ICOUP1(I) = ICOUP1(J+1)+ICPI
                        call COMP(I,LJ,ITYP,I,IT1,IT2)
                      end do
                      first3 = .false.
                    else
                      if (IT1 >= IT2) then
                        call LOOP5(KM,ISTOP,IT1,IT2)
                      else
                        call LOOP6(KM,ISTOP,IT1,IT2)
                      end if
                      if (ISTOP == 0) cycle loop_5
                    end if
                    KM = KM+1
                    if (KM == J) cycle loop_4
                  end do
                  if (.true.) exit loop_5
                end do loop_5
                exit loop_4
              else
                call PATH(KM,ISTOP,IT1,IT2)
                if (ISTOP == 0) then
                  skip = .true.
                  cycle loop_3
                end if
                KM = KM+1
                if (KM == K) cycle loop_3
              end if
            end do loop_4
            if (.true.) exit loop_3
          end do loop_3
        end if
      end if
      call LOOP5(KM,ISTOP,IT1,IT2)
      if (ISTOP == 0) then
        first1 = .true.
        cycle loop_2
      end if
      KM = KM+1
      if (KM == L) cycle loop_1
    end do loop_2
    if (.true.) exit loop_1
  end do loop_1
end do

return

end subroutine INT1
