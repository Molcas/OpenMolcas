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

subroutine INT7(I,K,L,IDIAG,BUFOUT,INDOUT,ICAD,IBUFL,KBUF,NTPB)
! I < L  I == K  J == L

use guga_global, only: COUP, IADD11, ICOUP, ICOUP1, ILIM, IJ, IRC, IV0, IVF0, IWAY, IX, J1, J2, JM, JM1, JNDX, LNP, Lu_11, MXVERT
use Constants, only: Half
use Definitions, only: wp, iwp, RtoI

implicit none
integer(kind=iwp), intent(in) :: I, K, L, IDIAG, ICAD(*), KBUF, NTPB
real(kind=wp), intent(inout) :: BUFOUT(*)
integer(kind=iwp), intent(inout) :: INDOUT(*), IBUFL(*)
integer(kind=iwp) :: IAD110, ICP, ICPP, ICQ, IDIV, IJJ, IN_, IPOS, ISTOP, ISU, ISUM, IT1, IT2, ITAIL, ITT, ITURN, ITYP, IVL, JND1, &
                     KBUF0, KBUF1, KBUF2, KM, LJ, LJM, LJS, NBN
logical(kind=iwp) :: first

IJJ = 0 ! dummy initialize
KBUF0 = RtoI*KBUF
KBUF1 = KBUF0+KBUF+1
KBUF2 = KBUF1+1
IDIV = RtoI
ITYP = 0
if (IDIAG == 1) IJJ = L*(L-1)/2+K
LJS = IJ(L+1)+1
LJM = IJ(L)
do ITT=1,ILIM
  IT1 = (ITT-1)*MXVERT
  IT2 = IT1
  do LJ=LJS,LJM
    ITURN = 0
    if (IDIAG == 1) ITURN = 1
    IWAY(L) = 1
    loop_1: do
      KM = L
      J2(KM+1) = LJ
      J1(KM+1) = LJ
      JM(KM) = IVF0+1
      JM1(KM) = IVF0+1
      if (ITURN == 0) then
        call LOOP7(KM,ISTOP,IT1,IT2)
      else if (ITURN == 1) then
        call LOOP8(KM,ISTOP,IT1,IT2)
      end if
      if (ISTOP /= 1) then
        if ((IDIAG == 1) .and. (J1(KM) /= J2(KM))) cycle loop_1
      else
        if (ITURN == 1) exit loop_1
        ITURN = 1
        IWAY(L) = 1
        cycle loop_1
      end if
      loop_2: do
        KM = KM-1
        IWAY(KM) = 1
        first = KM == I
        loop_3: do
          if (first) then
            KM = I
            if (ITURN == 0) then
              call LOOP14(KM,ISTOP,IT1,IT2)
            else if (ITURN == 1) then
              call LOOP18(KM,ISTOP,IT1,IT2)
            end if
            if (ISTOP == 1) then
              KM = KM+1
              if (KM == L) cycle loop_1
              first = .false.
              cycle loop_3
            end if
            if (abs(COUP(I)) < 1.0e-6_wp) cycle loop_3
            if (IDIAG /= 0) then
              if (ICOUP1(I) /= ICOUP(I)) cycle loop_3
            else if ((ITURN /= 0) .and. (IWAY(L) /= 5)) then
              if (ICOUP1(I) <= ICOUP(I)) cycle loop_3
            end if
            if (ITURN == 0) COUP(I) = Half*COUP(I)
            if (IDIAG /= 1) then
              call COMP(I,LJ,ITYP,I,IT1,IT2)
              cycle loop_3
            end if
            loop_4: do
              KM = KM-1
              if (KM /= 0) IWAY(KM) = 1
              do
                if (KM == 0) then
                  IVL = J2(1)
                  ITAIL = IX(IT1+LJ)
                  ISUM = IV0-IVL
                  ISU = 0
                  if (ISUM /= 0) ISU = IRC(ISUM)
                  do IN_=1,ITAIL
                    JND1 = JNDX(ISU+ICOUP(1)+IN_)
                    if (JND1 == 0) cycle
                    IPOS = (JND1-1)*LNP+IJJ
                    NBN = (IPOS-1)/NTPB+1
                    IBUFL(NBN) = IBUFL(NBN)+1
                    ICQ = ICAD(NBN)
                    ICP = ICQ/IDIV+IBUFL(NBN)
                    BUFOUT(ICP) = COUP(I)
                    ICPP = ICQ+KBUF0+IBUFL(NBN)
                    INDOUT(ICPP) = IPOS
                    if (IBUFL(NBN) < KBUF) cycle
                    INDOUT(ICQ+KBUF1) = KBUF
                    IAD110 = IADD11
                    call iDAFILE(Lu_11,1,INDOUT(ICQ+1),KBUF2,IADD11)
                    INDOUT(ICQ+KBUF2) = IAD110
                    IBUFL(NBN) = 0
                  end do
                  if (I == 1) cycle loop_3
                  KM = 1
                else
                  call PATH(KM,ISTOP,IT1,IT2)
                  if (ISTOP == 0) cycle loop_4
                  KM = KM+1
                  if (KM == I) cycle loop_3
                end if
              end do
              if (.true.) exit loop_4
            end do loop_4
            exit loop_3
          else
            JM(KM) = IVF0+1
            JM1(KM) = IVF0+1
            if (ITURN == 0) then
              call LOOP17(KM,ISTOP,IT1,IT2)
            else if (ITURN == 1) then
              call LOOP21(KM,ISTOP,IT1,IT2)
            end if
            if (ISTOP /= 1) then
              if ((IDIAG == 1) .and. (J1(KM) /= J2(KM))) cycle loop_3
              cycle loop_2
            end if
            KM = KM+1
            if (KM == L) cycle loop_1
          end if
        end do loop_3
        if (.true.) exit loop_2
      end do loop_2
      if (.true.) exit loop_1
    end do loop_1
  end do
end do

return

end subroutine INT7
