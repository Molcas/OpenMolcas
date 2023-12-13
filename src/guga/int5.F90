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

subroutine INT5(I,J,L)
! I < J < L   I == K

use guga_global, only: COUP, IJ, ILIM, IVF0, IWAY, J1, J2, JM, JM1, MXVERT
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: I, J, L
integer(kind=iwp) :: IFAI, ISTOP, IT1, IT2, ITT, ITURN, ITYP, KM, LJ, LJM, LJS
logical(kind=iwp) :: first1, first2, skip

ITYP = 0
LJS = IJ(L+1)+1
LJM = IJ(L)
do ITT=1,ILIM
  IT1 = (ITT-1)*MXVERT
  IT2 = IT1
  do LJ=LJS,LJM
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
            JM(KM) = IVF0+1
            JM1(KM) = IVF0+1
            if (ITURN == 0) then
              call LOOP10(KM,ISTOP,IT1,IT2)
            else if (ITURN == 1) then
              IFAI = 0
              call LOOP13(KM,ISTOP,IFAI,IT1,IT2)
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
              IWAY(KM) = 1
              first2 = KM == I
              do
                if (first2) then
                  do
                    KM = I
                    if (ITURN == 0) then
                      call LOOP14(KM,ISTOP,IT1,IT2)
                    else if (ITURN == 1) then
                      call LOOP22(KM,ISTOP,IT1,IT2)
                    end if
                    if (ISTOP == 1) exit
                    if (abs(COUP(I)) < 1.0e-6_wp) cycle
                    call COMP(I,LJ,ITYP,I,IT1,IT2)
                  end do
                  first2 = .false.
                else
                  JM(KM) = IVF0+1
                  JM1(KM) = IVF0+1
                  if (ITURN == 0) then
                    call LOOP17(KM,ISTOP,IT1,IT2)
                  else if (ITURN == 1) then
                    IFAI = 0
                    call LOOP23(KM,ISTOP,IFAI,IT1,IT2)
                  end if
                  if (ISTOP == 0) cycle loop_4
                end if
                KM = KM+1
                if (KM == J) cycle loop_3
              end do
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
end do

return

end subroutine INT5
