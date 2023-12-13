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

subroutine SS_HEAD_DBL_TAIL_ACT_G(LRA)

use gugaci_global, only: jb_sys, jml, jmr, jpad, jpadl, jpel, jper, just, jwl, jwr, line, lrs, lsm_inn, norb_dz, norb_frz, w0, &
                         w0_ss, w1, w1_ss
use Symmetry_Info, only: Mul
use Constants, only: Zero
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: lra
integer(kind=iwp) :: iwdl, iwdr, jmlr, list0, list0_1, list0_2, list0_3, list0_4, list1, list1_1, list1_2, list1_3, list1_4, lmi, &
                     lmij, lmj, lmk, lmki, lmkj, lr0, lri, lrj, lrk, ni
real(kind=wp) :: vlop0, vlop1, w0ss1, w0ss10, w0ss11, w0ss12, w0ss13, w0ss14, w0ss15, w0ss16, w0ss17, w0ss18, w0ss2, w0ss20, &
                 w0ss3, w0ss4, w0ss5, w0ss6, w0ss7, w0ss8, w0ss9, w1ss1, w1ss10, w1ss11, w1ss12, w1ss13, w1ss14, w1ss15, w1ss16, &
                 w1ss17, w1ss18, w1ss2, w1ss3, w1ss4, w1ss5, w1ss6, w1ss7, w1ss8, w1ss9, wl0, wl0_1, wl0_2, wl0_3, wl0_4, wl1, &
                 wl1_1, wl1_2, wl1_3, wl1_4, xwl1_1

w0ss1 = Zero
w1ss1 = Zero
w0ss2 = Zero
w1ss2 = Zero
w0ss3 = Zero
w1ss3 = Zero
w0ss4 = Zero
w1ss4 = Zero
w0ss5 = Zero
w1ss5 = Zero
w0ss6 = Zero
w1ss6 = Zero
w0ss7 = Zero
w1ss7 = Zero
w0ss8 = Zero
w1ss8 = Zero
w0ss9 = Zero
w1ss9 = Zero
w0ss10 = Zero
w1ss10 = Zero
w0ss11 = Zero
w1ss11 = Zero
w0ss12 = Zero
w1ss12 = Zero
w0ss13 = Zero
w1ss13 = Zero
w0ss14 = Zero
w1ss14 = Zero
w0ss15 = Zero
w1ss15 = Zero
w0ss16 = Zero
w1ss16 = Zero
w0ss18 = Zero
w1ss18 = Zero

! SS(1-1)  Ar(01)-Bl(32)-
! SS(1-2)  Ar(02)-Bl(31)-
! SS(1-3)  Ar(13)-Bl(20)-
! SS(1-4)  Ar(23)-Bl(10)-
! SS(1-5)  (22)-Ar(13)-Bl(31)-
! SS(1-6)  (11)-Ar(23)-Bl(32)-
! SS(1-7)  Ar(13)-C'(21)-Bl(32)-
! SS(1-8)  Ar(13)-C'(22)-Bl(31)-
! SS(1-9)  Ar(23)-C'(11)-Bl(32)-
! SS(1-10) Ar(23)-C'(12)-Bl(31)-
! SS(1-11) Ar(13)-Bl(31)-C"(22)-
! SS(1-12) Ar(13)-Bl(32)-C"(21)-
! SS(1-13) Ar(23)-Bl(31)-C"(12)-
! SS(1-14) Ar(23)-Bl(32)-C"(11)-
! SS(1-15) (22)-Drl(11)-
! SS(1-16) (11)-Drl(22)-
! SS(1-17) Drl(22)-C"(11)-
! SS(1-18) Drl(11)-C"(22)-
! SS(1-19) Drl(12)-C"(21)-
! SS(1-20) Drl(33)-C"(00)-
! SS(1-20) Drl(33)-C"(11)-C"(22)-
! SS(1-20) (11)Drl(33)-C"(22)-
! SS(1-20) (11)(22)Drl(33)-
! SS(1-20) Drl(33)-C"(22)-C"(11)-
! SS(1-20) (22)Drl(33)-C"(11)-
! SS(1-20) (22)(11)Drl(33)-
! ss(1-20) (21)(12)Drl(33)-

do LRI=NORB_FRZ+1,NORB_DZ-1
  LMI = LSM_INN(LRI)
  do LRJ=LRI+1,NORB_DZ
    LMJ = LSM_INN(LRJ)
    LMIJ = Mul(LMI,LMJ)
    JMLR = Mul(JML,JMR)
    if (LMIJ /= JMLR) cycle
    W0SS2 = W0_SS(2)
    W1SS2 = W1_SS(2)
    W0SS4 = W0_SS(4)
    W1SS4 = W1_SS(4)
    W0SS5 = W0_SS(5)
    W1SS5 = W1_SS(5)
    W0SS10 = -W0_SS(10)
    W1SS10 = -W1_SS(10)
    W0SS14 = W0_SS(14)
    W1SS14 = W1_SS(14)
    if (JB_SYS > 0) then
      W0SS1 = W0_SS(1)
      W1SS1 = W1_SS(1)
      W0SS3 = W0_SS(3)
      W1SS3 = W1_SS(3)
      W0SS6 = W0_SS(6)
      W1SS6 = W1_SS(6)
      W0SS7 = W0_SS(7)
      W1SS7 = W1_SS(7)
      W0SS8 = W0_SS(8)
      W1SS8 = W1_SS(8)
      W0SS9 = W0_SS(9)
      W1SS9 = W1_SS(9)
      W0SS11 = W0_SS(11)
      W1SS11 = W1_SS(11)
      W0SS12 = W0_SS(12)
      W1SS12 = W1_SS(12)
      W0SS13 = W0_SS(13)
      W1SS13 = W1_SS(13)
    end if
    NI = mod(LRJ-LRI,2)
    if (NI == 0) then
      W0SS2 = -W0SS2
      W1SS2 = -W1SS2
      W0SS4 = -W0SS4
      W1SS4 = -W1SS4
      W0SS5 = -W0SS5
      W1SS5 = -W1SS5
      W0SS10 = -W0SS10
      W1SS10 = -W1SS10
      W0SS14 = -W0SS14
      W1SS14 = -W1SS14
      if (JB_SYS > 0) then
        W0SS1 = -W0SS1
        W1SS1 = -W1SS1
        W0SS3 = -W0SS3
        W1SS3 = -W1SS3
        W0SS6 = -W0SS6
        W1SS6 = -W1SS6
        W0SS7 = -W0SS7
        W1SS7 = -W1SS7
        W0SS8 = -W0SS8
        W1SS8 = -W1SS8
        W0SS9 = -W0SS9
        W1SS9 = -W1SS9
        W0SS11 = -W0SS11
        W1SS11 = -W1SS11
        W0SS12 = -W0SS12
        W1SS12 = -W1SS12
        W0SS13 = -W0SS13
        W1SS13 = -W1SS13
      end if
    end if

    if ((JML == 1) .and. (LMIJ == JMR)) then
      IWDL = JUST(LRI,LRI)
      ! SS(1-1)   Ar(01)-Bl(32)-
      if (JB_SYS > 0) then
        IWDR = JUST(LRJ,LRI)
        VLOP0 = W0*W0SS1
        VLOP1 = W1*W1SS1
        if (LINE == 26) then    !LRI,LRJ,LRA
          call COMP_LOOP_G(5,LRI,LRJ,LRS,LRA,VLOP0,VLOP1,WL0,LIST0,WL1,LIST1)
        end if
        if (LINE == 28) then    !LRI,LRJ,LRS,LRA
          call COMP_LOOP_G(7,LRI,LRJ,LRS,LRA,VLOP0,VLOP1,WL0,LIST0,WL1,LIST1)
        end if
        if (LINE == 29) then    !LRI,LRJ,LRS,LRA
          call COMP_LOOP_G(6,LRI,LRJ,LRS,LRA,VLOP0,VLOP1,WL0,LIST0,WL1,LIST1)
        end if

        call PRODAB_2(3,JPEL,IWDL,IWDR,JWL,JWR,WL0,JPER,LIST0)
        if (LIST1 /= 0) call PRODAB_2(3,JPEL,IWDL,IWDR,JWL,JWR,WL1,JPER,LIST1)

      end if
      ! SS(1-2)  Ar(02)-Bl(31)-
      IWDR = JUST(LRI,LRJ)
      VLOP0 = W0*W0SS2
      VLOP1 = W1*W1SS2
      if (LINE == 26) then    !LRI,LRJ,LRA
        call COMP_LOOP_G(5,LRI,LRJ,LRS,LRA,VLOP0,VLOP1,WL0,LIST0,WL1,LIST1)
      end if
      if (LINE == 28) then    !LRI,LRJ,LRS,LRA
        call COMP_LOOP_G(7,LRI,LRJ,LRS,LRA,VLOP0,VLOP1,WL0,LIST0,WL1,LIST1)
      end if
      if (LINE == 29) then    !LRI,LRJ,LRS,LRA
        call COMP_LOOP_G(6,LRI,LRJ,LRS,LRA,VLOP0,VLOP1,WL0,LIST0,WL1,LIST1)
      end if
      call PRODAB_2(3,JPEL,IWDL,IWDR,JWL,JWR,WL0,JPER,LIST0)
      if (LIST1 /= 0) call PRODAB_2(3,JPEL,IWDL,IWDR,JWL,JWR,WL1,JPER,LIST1)

    end if

    if ((JMR == 1) .and. (LMIJ == JML)) then
      IWDR = JUST(LRJ,LRJ)
      ! SS(1-3)  Ar(13)-Bl(20)-
      if (JB_SYS > 0) then
        IWDL = JUST(LRJ,LRI)
        VLOP0 = W0*W0SS3
        VLOP1 = W1*W1SS3
        if (LINE == 26) then    !LRI,LRJ,LRA
          call COMP_LOOP_G(5,LRI,LRJ,LRS,LRA,VLOP0,VLOP1,WL0,LIST0,WL1,LIST1)
        end if
        if (LINE == 28) then    !LRI,LRJ,LRS,LRA
          call COMP_LOOP_G(7,LRI,LRJ,LRS,LRA,VLOP0,VLOP1,WL0,LIST0,WL1,LIST1)
        end if
        if (LINE == 29) then    !LRI,LRJ,LRS,LRA
          call COMP_LOOP_G(6,LRI,LRJ,LRS,LRA,VLOP0,VLOP1,WL0,LIST0,WL1,LIST1)
        end if
        call PRODAB_2(3,JPEL,IWDL,IWDR,JWL,JWR,WL0,JPER,LIST0)
        if (LIST1 /= 0) call PRODAB_2(3,JPEL,IWDL,IWDR,JWL,JWR,WL1,JPER,LIST1)

      end if
      ! SS(1-4)  Ar(23)-Bl(10)-        ACT -C"-                  ! IPRAD
      IWDL = JUST(LRI,LRJ)
      VLOP0 = W0*W0SS4
      VLOP1 = W1*W1SS4
      if (LINE == 26) then    !LRI,LRJ,LRA
        call COMP_LOOP_G(5,LRI,LRJ,LRS,LRA,VLOP0,VLOP1,WL0,LIST0,WL1,LIST1)
      end if
      if (LINE == 28) then    !LRI,LRJ,LRS,LRA
        call COMP_LOOP_G(7,LRI,LRJ,LRS,LRA,VLOP0,VLOP1,WL0,LIST0,WL1,LIST1)
      end if
      if (LINE == 29) then    !LRI,LRJ,LRS,LRA
        call COMP_LOOP_G(6,LRI,LRJ,LRS,LRA,VLOP0,VLOP1,WL0,LIST0,WL1,LIST1)
      end if
      call PRODAB_2(3,JPEL,IWDL,IWDR,JWL,JWR,WL0,JPER,LIST0)
      if (LIST1 /= 0) call PRODAB_2(3,JPEL,IWDL,IWDR,JWL,JWR,WL1,JPER,LIST1)
    end if

    VLOP0 = W0*W0SS5
    VLOP1 = W1*W1SS5
    if (LINE == 26) then    !LRI,LRJ,LRA
      call COMP_LOOP_G(5,LRI,LRJ,LRS,LRA,VLOP0,VLOP1,WL0,LIST0,WL1,LIST1)
    end if
    if (LINE == 28) then    !LRI,LRJ,LRS,LRA
      call COMP_LOOP_G(7,LRI,LRJ,LRS,LRA,VLOP0,VLOP1,WL0,LIST0,WL1,LIST1)
    end if
    if (LINE == 29) then    !LRI,LRJ,LRS,LRA
      call COMP_LOOP_G(6,LRI,LRJ,LRS,LRA,VLOP0,VLOP1,WL0,LIST0,WL1,LIST1)
    end if
    if (JB_SYS > 0) then
      VLOP0 = W0*W0SS6
      VLOP1 = W1*W1SS6
      if (LINE == 26) then    !LRI,LRJ,LRA
        call COMP_LOOP_G(5,LRI,LRJ,LRS,LRA,VLOP0,VLOP1,WL0_1,LIST0_1,WL1_1,LIST1_1)
      end if
      if (LINE == 28) then    !LRI,LRJ,LRS,LRA
        call COMP_LOOP_G(7,LRI,LRJ,LRS,LRA,VLOP0,VLOP1,WL0_1,LIST0_1,WL1_1,LIST1_1)
      end if
      if (LINE == 29) then    !LRI,LRJ,LRS,LRA
        call COMP_LOOP_G(6,LRI,LRJ,LRS,LRA,VLOP0,VLOP1,WL0_1,LIST0_1,WL1_1,LIST1_1)
      end if
    end if
    do LRK=NORB_FRZ+1,LRI-1
      LMK = LSM_INN(LRK)
      LMKI = Mul(LMK,LMI)
      LMKJ = Mul(LMK,LMJ)
      if ((LMKI == JML) .and. (LMKJ == JMR)) then
        ! SS(1-5)  (22)-Ar(13)-Bl(31)-
        IWDL = JUST(LRK,LRI)
        IWDR = JUST(LRK,LRJ)
        call PRODAB_2(3,JPEL,IWDL,IWDR,JWL,JWR,WL0,JPER,LIST0)
        if (LIST1 /= 0) call PRODAB_2(3,JPEL,IWDL,IWDR,JWL,JWR,WL1,JPER,LIST1)

        if (JB_SYS > 0) then
          ! SS(1-6)  (11)-Ar(23)-Bl(32)-
          IWDL = JUST(LRI,LRK)
          IWDR = JUST(LRJ,LRK)
          call PRODAB_2(3,JPEL,IWDL,IWDR,JWL,JWR,WL0_1,JPER,LIST0_1)
          if (LIST1_1 /= 0) call PRODAB_2(3,JPEL,IWDL,IWDR,JWL,JWR,WL1_1,JPER,LIST1_1)

        end if
      end if
    end do

    if (JB_SYS > 0) then
      VLOP0 = W0*W0SS7
      VLOP1 = W1*W1SS7
      if (LINE == 26) then    !LRI,LRJ,LRA
        call COMP_LOOP_G(5,LRI,LRJ,LRS,LRA,VLOP0,VLOP1,WL0_1,LIST0_1,WL1_1,LIST1_1)
      end if
      if (LINE == 28) then    !LRI,LRJ,LRS,LRA
        call COMP_LOOP_G(7,LRI,LRJ,LRS,LRA,VLOP0,VLOP1,WL0_1,LIST0_1,WL1_1,LIST1_1)
      end if
      if (LINE == 29) then    !LRI,LRJ,LRS,LRA
        call COMP_LOOP_G(6,LRI,LRJ,LRS,LRA,VLOP0,VLOP1,WL0_1,LIST0_1,WL1_1,LIST1_1)
      end if
      VLOP0 = W0*W0SS8
      VLOP1 = W1*W1SS8
      if (LINE == 26) then    !LRI,LRJ,LRA
        call COMP_LOOP_G(5,LRI,LRJ,LRS,LRA,VLOP0,VLOP1,WL0_2,LIST0_2,WL1_2,LIST1_2)
      end if
      if (LINE == 28) then    !LRI,LRJ,LRS,LRA
        call COMP_LOOP_G(7,LRI,LRJ,LRS,LRA,VLOP0,VLOP1,WL0_2,LIST0_2,WL1_2,LIST1_2)
      end if
      if (LINE == 29) then    !LRI,LRJ,LRS,LRA
        call COMP_LOOP_G(6,LRI,LRJ,LRS,LRA,VLOP0,VLOP1,WL0_2,LIST0_2,WL1_2,LIST1_2)
      end if
      VLOP0 = W0*W0SS9
      VLOP1 = W1*W1SS9
      if (LINE == 26) then    !LRI,LRJ,LRA
        call COMP_LOOP_G(5,LRI,LRJ,LRS,LRA,VLOP0,VLOP1,WL0_3,LIST0_3,WL1_3,LIST1_3)
      end if
      if (LINE == 28) then    !LRI,LRJ,LRS,LRA
        call COMP_LOOP_G(7,LRI,LRJ,LRS,LRA,VLOP0,VLOP1,WL0_3,LIST0_3,WL1_3,LIST1_3)
      end if
      if (LINE == 29) then    !LRI,LRJ,LRS,LRA
        call COMP_LOOP_G(6,LRI,LRJ,LRS,LRA,VLOP0,VLOP1,WL0_3,LIST0_3,WL1_3,LIST1_3)
      end if
    end if
    VLOP0 = W0*W0SS10
    VLOP1 = W1*W1SS10
    if (LINE == 26) then    !LRI,LRJ,LRA
      call COMP_LOOP_G(5,LRI,LRJ,LRS,LRA,VLOP0,VLOP1,WL0_4,LIST0_4,WL1_4,LIST1_4)
    end if
    if (LINE == 28) then    !LRI,LRJ,LRS,LRA
      call COMP_LOOP_G(7,LRI,LRJ,LRS,LRA,VLOP0,VLOP1,WL0_4,LIST0_4,WL1_4,LIST1_4)
    end if
    if (LINE == 29) then    !LRI,LRJ,LRS,LRA
      call COMP_LOOP_G(6,LRI,LRJ,LRS,LRA,VLOP0,VLOP1,WL0_4,LIST0_4,WL1_4,LIST1_4)
    end if

    do LRK=LRI+1,LRJ-1
      LMK = LSM_INN(LRK)
      LMKI = Mul(LMK,LMI)
      LMKJ = Mul(LMK,LMJ)
      if ((LMKI == JML) .and. (LMKJ == JMR)) then
        if (JB_SYS > 0) then
          ! SS(1-7)  Ar(13)-C'(21)-Bl(32)-
          IWDL = JUST(LRK,LRI)
          IWDR = JUST(LRJ,LRK)
          !call PRODAB(3,JPEL,IWDL,IWDR,JWL,JWR,-WL1,JPER)
          call PRODAB_2(3,JPEL,IWDL,IWDR,JWL,JWR,-WL0_1,JPER,LIST0_1)
          if (LIST1_1 /= 0) call PRODAB_2(3,JPEL,IWDL,IWDR,JWL,JWR,-WL1_1,JPER,LIST1_1)

          ! SS(1-8)  Ar(13)-C'(22)-Bl(31)-
          IWDL = JUST(LRK,LRI)
          IWDR = JUST(LRK,LRJ)
          !call PRODAB(3,JPEL,IWDL,IWDR,JWL,JWR,-WL2,JPER)
          call PRODAB_2(3,JPEL,IWDL,IWDR,JWL,JWR,-WL0_2,JPER,LIST0_2)
          if (LIST1_2 /= 0) call PRODAB_2(3,JPEL,IWDL,IWDR,JWL,JWR,-WL1_2,JPER,LIST1_2)

          ! SS(1-9)  Ar(23)-C'(11)-Bl(32)-
          IWDL = JUST(LRI,LRK)
          IWDR = JUST(LRJ,LRK)
          !WL = (VLOP0-VLOP1)*vint_ci(LIST)-2*VLOP0*vint_ci(LIST+1)
          !call PRODAB(3,JPEL,IWDL,IWDR,JWL,JWR,-WL3,JPER)
          call PRODAB_2(3,JPEL,IWDL,IWDR,JWL,JWR,-WL0_3,JPER,LIST0_3)
          if (LIST1_3 /= 0) call PRODAB_2(3,JPEL,IWDL,IWDR,JWL,JWR,-WL1_3,JPER,LIST1_3)

        end if
        ! SS(1-10) Ar(23)-C'(12)-Bl(31)-
        IWDL = JUST(LRI,LRK)
        IWDR = JUST(LRK,LRJ)
        !call PRODAB(3,JPEL,IWDL,IWDR,JWL,JWR,WL4,JPER)
        call PRODAB_2(3,JPEL,IWDL,IWDR,JWL,JWR,WL0_4,JPER,LIST0_4)
        if (LIST1_4 /= 0) call PRODAB_2(3,JPEL,IWDL,IWDR,JWL,JWR,WL1_4,JPER,LIST1_4)

      end if
    end do

    if (JB_SYS > 0) then
      VLOP0 = W0*W0SS11
      VLOP1 = W1*W1SS11
      if (LINE == 26) then    !LRI,LRJ,LRA
        call COMP_LOOP_G(5,LRI,LRJ,LRS,LRA,VLOP0,VLOP1,WL0_1,LIST0_1,xWL1_1,LIST1_1)
      end if
      if (LINE == 28) then    !LRI,LRJ,LRS,LRA
        call COMP_LOOP_G(7,LRI,LRJ,LRS,LRA,VLOP0,VLOP1,WL0_1,LIST0_1,WL1_1,LIST1_1)
      end if
      if (LINE == 29) then    !LRI,LRJ,LRS,LRA
        call COMP_LOOP_G(6,LRI,LRJ,LRS,LRA,VLOP0,VLOP1,WL0_1,LIST0_1,WL1_1,LIST1_1)
      end if
      VLOP0 = W0*W0SS12
      VLOP1 = W1*W1SS12
      if (LINE == 26) then    !LRI,LRJ,LRA
        call COMP_LOOP_G(5,LRI,LRJ,LRS,LRA,VLOP0,VLOP1,WL0_2,LIST0_2,WL1_2,LIST1_2)
      end if
      if (LINE == 28) then    !LRI,LRJ,LRS,LRA
        call COMP_LOOP_G(7,LRI,LRJ,LRS,LRA,VLOP0,VLOP1,WL0_2,LIST0_2,WL1_2,LIST1_2)
      end if
      if (LINE == 29) then    !LRI,LRJ,LRS,LRA
        call COMP_LOOP_G(6,LRI,LRJ,LRS,LRA,VLOP0,VLOP1,WL0_2,LIST0_2,WL1_2,LIST1_2)
      end if
      VLOP0 = W0*W0SS13
      VLOP1 = W1*W1SS13
      if (LINE == 26) then    !LRI,LRJ,LRA
        call COMP_LOOP_G(5,LRI,LRJ,LRS,LRA,VLOP0,VLOP1,WL0_3,LIST0_3,WL1_3,LIST1_3)
      end if
      if (LINE == 28) then    !LRI,LRJ,LRS,LRA
        call COMP_LOOP_G(7,LRI,LRJ,LRS,LRA,VLOP0,VLOP1,WL0_3,LIST0_3,WL1_3,LIST1_3)
      end if
      if (LINE == 29) then    !LRI,LRJ,LRS,LRA
        call COMP_LOOP_G(6,LRI,LRJ,LRS,LRA,VLOP0,VLOP1,WL0_3,LIST0_3,WL1_3,LIST1_3)
      end if
    end if
    VLOP0 = W0*W0SS14
    VLOP1 = W1*W1SS14
    if (LINE == 26) then    !LRI,LRJ,LRA
      call COMP_LOOP_G(5,LRI,LRJ,LRS,LRA,VLOP0,VLOP1,WL0_4,LIST0_4,WL1_4,LIST1_4)
    end if
    if (LINE == 28) then    !LRI,LRJ,LRS,LRA
      call COMP_LOOP_G(7,LRI,LRJ,LRS,LRA,VLOP0,VLOP1,WL0_4,LIST0_4,WL1_4,LIST1_4)
    end if
    if (LINE == 29) then    !LRI,LRJ,LRS,LRA
      call COMP_LOOP_G(6,LRI,LRJ,LRS,LRA,VLOP0,VLOP1,WL0_4,LIST0_4,WL1_4,LIST1_4)
    end if

    do LRK=LRJ+1,NORB_DZ
      LMK = LSM_INN(LRK)
      LMKI = Mul(LMK,LMI)
      LMKJ = Mul(LMK,LMJ)
      if ((LMKI == JML) .and. (LMKJ == JMR)) then
        if (JB_SYS > 0) then
          ! SS(1-11) Ar(13)-Bl(31)-C"(22)-
          IWDL = JUST(LRK,LRI)
          IWDR = JUST(LRK,LRJ)
          !call PRODAB(3,JPEL,IWDL,IWDR,JWL,JWR,WL1,JPER)
          call PRODAB_2(3,JPEL,IWDL,IWDR,JWL,JWR,WL0_1,JPER,LIST0_1)
          if (LIST1_1 /= 0) call PRODAB_2(3,JPEL,IWDL,IWDR,JWL,JWR,WL1_1,JPER,LIST1_1)

          ! SS(1-12) Ar(13)-Bl(32)-C"(21)-
          IWDL = JUST(LRK,LRI)
          IWDR = JUST(LRJ,LRK)
          !call PRODAB(3,JPEL,IWDL,IWDR,JWL,JWR,WL2,JPER)
          call PRODAB_2(3,JPEL,IWDL,IWDR,JWL,JWR,WL0_2,JPER,LIST0_2)
          if (LIST1_2 /= 0) call PRODAB_2(3,JPEL,IWDL,IWDR,JWL,JWR,WL1_2,JPER,LIST1_2)

          ! SS(1-13) Ar(23)-Bl(31)-C"(12)-
          IWDL = JUST(LRI,LRK)
          IWDR = JUST(LRK,LRJ)
          !call PRODAB(3,JPEL,IWDL,IWDR,JWL,JWR,WL3,JPER)
          call PRODAB_2(3,JPEL,IWDL,IWDR,JWL,JWR,WL0_3,JPER,LIST0_3)
          if (LIST1_3 /= 0) call PRODAB_2(3,JPEL,IWDL,IWDR,JWL,JWR,WL1_3,JPER,LIST1_3)

        end if
        ! SS(1-14) Ar(23)-Bl(32)-C"(11)- ACT -C"-
        IWDL = JUST(LRI,LRK)
        IWDR = JUST(LRJ,LRK)
        !call PRODAB(3,JPEL,IWDL,IWDR,JWL,JWR,WL4,JPER)
        call PRODAB_2(3,JPEL,IWDL,IWDR,JWL,JWR,WL0_4,JPER,LIST0_4)
        if (LIST1_4 /= 0) call PRODAB_2(3,JPEL,IWDL,IWDR,JWL,JWR,WL1_4,JPER,LIST1_4)

      end if
    end do
  end do
end do

if (JPAD /= JPADL) return

!if ((JB_SYS > 0) .or. (JWL >= JWR)) then
if (JB_SYS > 0) then
  if ((JPAD > 17) .and. (JPAD < 25)) then
    do LRI=NORB_FRZ+1,NORB_DZ
      LMI = LSM_INN(LRI)
      do LRJ=LRI+1,NORB_DZ
        LMJ = LSM_INN(LRJ)
        LMIJ = Mul(LMI,LMJ)
        if ((LMIJ /= JML) .or. (JML /= JMR)) cycle
        IWDL = JUST(LRJ,LRI)
        IWDR = JUST(LRI,LRJ)
        ! SS(1-19) Drl(12)-C"(21)-
        VLOP0 = W0*W0_SS(19)
        VLOP1 = W1*W1_SS(19)
        if (LINE == 26) then    !LRI,LRA
          call COMP_LOOP_G(9,LRI,0,0,LRA,VLOP0,VLOP1,WL0,LIST0,WL1,LIST1)
        end if
        if (LINE == 28) then    !LRI,LRS,LRA
          call COMP_LOOP_G(12,LRI,0,LRS,LRA,VLOP0,VLOP1,WL0,LIST0,WL1,LIST1)
        end if
        if (LINE == 29) then    !LRI,LRS,LRA
          call COMP_LOOP_G(11,LRI,0,LRS,LRA,VLOP0,VLOP1,WL0,LIST0,WL1,LIST1)
        end if
        !WL = (VLOP0-VLOP1)*VOINT(LRI,LRB)
        !call PRODAB(3,JPEL,IWDL,IWDR,JWL,JWR,WL,JPER)
        call PRODAB_2(3,JPEL,IWDL,IWDR,JWL,JWR,WL0,JPER,LIST0)
        if (LIST1 /= 0) call PRODAB_2(3,JPEL,IWDL,IWDR,JWL,JWR,WL1,JPER,LIST1)

      end do
    end do
  end if
end if
if (JWL >= JWR) return

W0SS15 = W0_SS(15)
W1SS15 = W1_SS(15)
W0SS17 = W0_SS(17)
W1SS17 = W1_SS(17)
W0SS20 = W0_SS(20)
if (JB_SYS > 0) then
  W0SS16 = W0_SS(16)
  W1SS16 = W1_SS(16)
  W0SS18 = W0_SS(18)
  W1SS18 = W1_SS(18)
end if

if ((JML == 1) .and. (JMR == 1)) then
  ! SS(1-20) Drl(33)-C"(00)-                         ! IPL(R)AD=1 or =NS
  do LR0=NORB_FRZ+1,NORB_DZ
    IWDL = JUST(LR0,LR0)
    IWDR = IWDL
    VLOP0 = W0*W0_SS(20)
    VLOP1 = Zero
    !WL = Zero
    do LRK=1,NORB_DZ
      if (LRK == LR0) cycle
      if (LINE == 26) then    !LRK,LRA
        call COMP_LOOP_G(9,LRK,0,0,LRA,VLOP0,VLOP1,WL0_1,LIST0_1,WL1_1,LIST1_1)
      end if
      if (LINE == 28) then    !LRK,LRS,LRA    LRS,LRA,LRK
        call COMP_LOOP_G(12,LRK,0,LRS,LRA,VLOP0,VLOP1,WL0_1,LIST0_1,WL1_1,LIST1_1)
      end if
      if (LINE == 29) then    !LRK,LRS,LRA
        call COMP_LOOP_G(11,LRK,0,LRS,LRA,VLOP0,VLOP1,WL0_1,LIST0_1,WL1_1,LIST1_1)
      end if
      !WL = WL+WL1
      call PRODAB_2(3,JPEL,IWDL,IWDR,JWL,JWR,WL0_1,JPER,LIST0_1)
      if (LIST1_1 /= 0) call PRODAB_2(3,JPEL,IWDL,IWDR,JWL,JWR,WL1_1,JPER,LIST1_1)

    end do
    !call PRODAB(3,JPEL,IWDL,IWDR,JWL,JWR,WL,JPER)
  end do
end if
do LRI=NORB_FRZ+1,NORB_DZ
  LMI = LSM_INN(LRI)
  do LRJ=LRI+1,NORB_DZ
    LMJ = LSM_INN(LRJ)
    LMIJ = Mul(LMI,LMJ)
    if ((LMIJ /= JML) .or. (JML /= JMR)) cycle
    if ((JWL >= JWR) .and. (JB_SYS == 0)) cycle
    !WL = Zero
    IWDL = JUST(LRI,LRJ)
    IWDR = IWDL
    ! SS(1-15) (22)-Drl(11)-
    VLOP0 = W0*W0SS15
    VLOP1 = W1*W1SS15
    !WL = Zero
    if (LINE == 26) then    !LRJ,LRA
      call COMP_LOOP_G(9,LRJ,0,0,LRA,VLOP0,VLOP1,WL0_1,LIST0_1,WL1_1,LIST1_1)
    end if
    if (LINE == 28) then    !LRJ,LRS,LRA
      call COMP_LOOP_G(12,LRJ,0,LRS,LRA,VLOP0,VLOP1,WL0_1,LIST0_1,WL1_1,LIST1_1)
    end if
    if (LINE == 29) then    !LRJ,LRS,LRA
      call COMP_LOOP_G(11,LRJ,0,LRS,LRA,VLOP0,VLOP1,WL0_1,LIST0_1,WL1_1,LIST1_1)
    end if
    !WL = WL+WL1
    call PRODAB_2(3,JPEL,IWDL,IWDR,JWL,JWR,WL0_1,JPER,LIST0_1)
    if (LIST1_1 /= 0) call PRODAB_2(3,JPEL,IWDL,IWDR,JWL,JWR,WL1_1,JPER,LIST1_1)

    !WL = WL+(VLOP0-VLOP1)*VOINT(LRJ,LRB)
    ! SS(1-17) Drl(22)-C"(11)-
    VLOP0 = W0*W0SS17
    VLOP1 = W1*W1SS17
    if (LINE == 26) then    !LRI,LRA
      call COMP_LOOP_G(9,LRI,0,0,LRA,VLOP0,VLOP1,WL0_1,LIST0_1,WL1_1,LIST1_1)
    end if
    if (LINE == 28) then    !LRI,LRS,LRA
      call COMP_LOOP_G(12,LRI,0,LRS,LRA,VLOP0,VLOP1,WL0_1,LIST0_1,WL1_1,LIST1_1)
    end if
    if (LINE == 29) then    !LRI,LRS,LRA
      call COMP_LOOP_G(11,LRI,0,LRS,LRA,VLOP0,VLOP1,WL0_1,LIST0_1,WL1_1,LIST1_1)
    end if
    !WL = WL+WL1
    call PRODAB_2(3,JPEL,IWDL,IWDR,JWL,JWR,WL0_1,JPER,LIST0_1)
    if (LIST1_1 /= 0) call PRODAB_2(3,JPEL,IWDL,IWDR,JWL,JWR,WL1_1,JPER,LIST1_1)

    !WL = WL+(VLOP0-VLOP1)*VOINT(LRI,LRB)
    ! SS(1-20) (22)(11)Drl(33)-
    ! SS(1-20) (22)Drl(33)-C"(11)-
    ! SS(1-20) Drl(33)-C"(22)-C"(11)-
    VLOP0 = W0*W0SS20
    VLOP1 = Zero
    do LRK=1,NORB_DZ
      if (LRK == LRI) cycle
      if (LRK == LRJ) cycle
      if (LINE == 26) then    !LRK,LRA
        call COMP_LOOP_G(9,LRK,0,0,LRA,VLOP0,VLOP1,WL0_1,LIST0_1,WL1_1,LIST1_1)
      end if
      if (LINE == 28) then    !LRK,LRS,LRA
        call COMP_LOOP_G(12,LRK,0,LRS,LRA,VLOP0,VLOP1,WL0_1,LIST0_1,WL1_1,LIST1_1)
      end if
      if (LINE == 29) then    !LRK,LRS,LRA
        call COMP_LOOP_G(11,LRK,0,LRS,LRA,VLOP0,VLOP1,WL0_1,LIST0_1,WL1_1,LIST1_1)
      end if
      !WL = WL+WL1
      call PRODAB_2(3,JPEL,IWDL,IWDR,JWL,JWR,WL0_1,JPER,LIST0_1)
      if (LIST1_1 /= 0) call PRODAB_2(3,JPEL,IWDL,IWDR,JWL,JWR,WL1_1,JPER,LIST1_1)

      !WL = WL+VLOP0*VOINT(LRK,LRB)
    end do
    !call PRODAB(3,JPEL,IWDL,IWDR,JWL,JWR,WL,JPER)
    if (JB_SYS > 0) then
      IWDL = JUST(LRJ,LRI)
      IWDR = IWDL
      ! SS(1-16) (11)-Drl(22)-
      VLOP0 = W0*W0SS16
      VLOP1 = W1*W1SS16
      !WL = Zero
      if (LINE == 26) then    !LRJ,LRA
        call COMP_LOOP_G(9,LRJ,0,0,LRA,VLOP0,VLOP1,WL0_1,LIST0_1,WL1_1,LIST1_1)
      end if
      if (LINE == 28) then    !LRJ,LRS,LRA
        call COMP_LOOP_G(12,LRJ,0,LRS,LRA,VLOP0,VLOP1,WL0_1,LIST0_1,WL1_1,LIST1_1)
      end if
      if (LINE == 29) then    !LRJ,LRS,LRA
        call COMP_LOOP_G(11,LRJ,0,LRS,LRA,VLOP0,VLOP1,WL0_1,LIST0_1,WL1_1,LIST1_1)
      end if
      !WL = WL+WL1
      call PRODAB_2(3,JPEL,IWDL,IWDR,JWL,JWR,WL0_1,JPER,LIST0_1)
      if (LIST1_1 /= 0) call PRODAB_2(3,JPEL,IWDL,IWDR,JWL,JWR,WL1_1,JPER,LIST1_1)

      !WL = (VLOP0-VLOP1)*VOINT(LRJ,LRB)
      ! SS(1-18) Drl(11)-C"(22)-
      VLOP0 = W0*W0SS18
      VLOP1 = W1*W1SS18
      if (LINE == 26) then    !LRI,LRA
        call COMP_LOOP_G(9,LRI,0,0,LRA,VLOP0,VLOP1,WL0_1,LIST0_1,WL1_1,LIST1_1)
      end if
      if (LINE == 28) then    !LRI,LRS,LRA
        call COMP_LOOP_G(12,LRI,0,LRS,LRA,VLOP0,VLOP1,WL0_1,LIST0_1,WL1_1,LIST1_1)
      end if
      if (LINE == 29) then    !LRI,LRS,LRA
        call COMP_LOOP_G(11,LRI,0,LRS,LRA,VLOP0,VLOP1,WL0_1,LIST0_1,WL1_1,LIST1_1)
      end if
      !WL = WL+WL1
      call PRODAB_2(3,JPEL,IWDL,IWDR,JWL,JWR,WL0_1,JPER,LIST0_1)
      if (LIST1_1 /= 0) call PRODAB_2(3,JPEL,IWDL,IWDR,JWL,JWR,WL1_1,JPER,LIST1_1)

      !WL = WL+(VLOP0-VLOP1)*VOINT(LRI,LRB)
      ! SS(1-20) Drl(33)-C"(11)-C"(22)-
      ! SS(1-20) (11)Drl(33)-C"(22)-
      ! SS(1-20) (11)(22)Drl(33)-
      VLOP0 = W0*W0SS20
      VLOP1 = Zero
      do LRK=1,NORB_DZ
        if (LRK == LRI) cycle
        if (LRK == LRJ) cycle
        if (LINE == 26) then    !LRK,LRA
          call COMP_LOOP_G(9,LRK,0,0,LRA,VLOP0,VLOP1,WL0_1,LIST0_1,WL1_1,LIST1_1)
        end if
        if (LINE == 28) then    !LRK,LRS,LRA
          call COMP_LOOP_G(12,LRK,0,LRS,LRA,VLOP0,VLOP1,WL0_1,LIST0_1,WL1_1,LIST1_1)
        end if
        if (LINE == 29) then    !LRK,LRS,LRA
          call COMP_LOOP_G(11,LRK,0,LRS,LRA,VLOP0,VLOP1,WL0_1,LIST0_1,WL1_1,LIST1_1)
        end if
        !WL = WL+WL1
        call PRODAB_2(3,JPEL,IWDL,IWDR,JWL,JWR,WL0_1,JPER,LIST0_1)
        if (LIST1_1 /= 0) call PRODAB_2(3,JPEL,IWDL,IWDR,JWL,JWR,WL1_1,JPER,LIST1_1)

        !WL = WL+VLOP0*VOINT(LRK,LRB)
      end do
      !call PRODAB(3,JPEL,IWDL,IWDR,JWL,JWR,WL,JPER)
    end if
  end do
end do

return

end subroutine SS_HEAD_DBL_TAIL_ACT_G

subroutine ST_HEAD_DBL_TAIL_ACT_G(LRA)

use gugaci_global, only: jb_sys, jml, jmr, jpel, jper, just, jwl, jwr, line, lrs, lsm_inn, norb_dz, norb_frz, w1, w1_st
use Symmetry_Info, only: Mul
use Constants, only: Zero
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: lra
integer(kind=iwp) :: iwds, iwdt, jmlr, list0, list0_1, list1, list1_1, lmi, lmij, lmj, lmk, lri, lrj, lrk, ni
real(kind=wp) :: vlop0, vlop1, w1st1, w1st2, w1st3, w1st4, wl0, wl0_1, wl1, wl1_1

! ST(2-1) Ar(02)-Bl(32)-
! ST(2-2) (22)Ar(13)-Bl(32)-
! ST(2-4) Ar(23)-C'(12)-Bl(32)-
! ST(2-4) Ar(23)-Bl(32)-C'(12)-
! ST(2-5) (22)Drl(12)-
! ST(2-6) Drl(22)-C"(12)-
do LRI=NORB_FRZ+1,NORB_DZ
  LMI = LSM_INN(LRI)
  do LRJ=LRI+1,NORB_DZ
    LMJ = LSM_INN(LRJ)
    LMIJ = Mul(LMI,LMJ)
    LMIJ = Mul(LMIJ,1)
    if ((JML == JMR) .and. (LMIJ == JML)) then
      IWDS = JUST(LRI,LRJ)
      IWDT = IWDS
      ! ST(2-5) (22)Drl(12)-
      VLOP1 = W1*W1_ST(5)             !D2-5
      VLOP0 = Zero
      if (LINE == 26) then    !LRJ,LRA
        call COMP_LOOP_G(9,LRJ,0,0,LRA,VLOP0,VLOP1,WL0,LIST0,WL1,LIST1)
      end if
      if (LINE == 28) then    !LRJ,LRS,LRA
        call COMP_LOOP_G(12,LRJ,0,LRS,LRA,VLOP0,VLOP1,WL0,LIST0,WL1,LIST1)
      end if
      if (LINE == 29) then    !LRJ,LRS,LRA
        call COMP_LOOP_G(11,LRJ,0,LRS,LRA,VLOP0,VLOP1,WL0,LIST0,WL1,LIST1)
      end if
      !call PRODAB_1(3,JPEL,IWDS,IWDT,JWL,JWR,WL2,JPER,LIST2,LIST3)
      call PRODAB_2(3,JPEL,IWDS,IWDT,JWL,JWR,WL0,JPER,LIST0)
      if (LIST1 /= 0) call PRODAB_2(3,JPEL,IWDS,IWDT,JWL,JWR,WL1,JPER,LIST1)

      ! ST(2-6) Drl(22)-C"(12)-
      VLOP1 = W1*W1_ST(6)             !D2-6
      VLOP0 = Zero
      if (LINE == 26) then    !LRI,LRA
        call COMP_LOOP_G(9,LRI,0,0,LRA,VLOP0,VLOP1,WL0_1,LIST0_1,WL1_1,LIST1_1)
      end if
      if (LINE == 28) then    !LRI,LRS,LRA
        call COMP_LOOP_G(12,LRI,0,LRS,LRA,VLOP0,VLOP1,WL0_1,LIST0_1,WL1_1,LIST1_1)
      end if
      if (LINE == 29) then    !LRI,LRS,LRA
        call COMP_LOOP_G(11,LRI,0,LRS,LRA,VLOP0,VLOP1,WL0_1,LIST0_1,WL1_1,LIST1_1)
      end if

      !call PRODAB_1(3,JPEL,IWDS,IWDT,JWL,JWR,WL2_1,JPER,LIST2_1,LIST3_1)
      call PRODAB_2(3,JPEL,IWDS,IWDT,JWL,JWR,WL0_1,JPER,LIST0_1)
      if (LIST1_1 /= 0) call PRODAB_2(3,JPEL,IWDS,IWDT,JWL,JWR,WL1_1,JPER,LIST1_1)

      !WL = WL+WL1
      !LIST = LIST3(LRA,LRB,LRI)
      !WL = WL-VLOP1*vint_ci(LIST)    !4.3 Vlop0=0        !!!!!
      !call PRODAB(3,JPEL,IWDS,IWDT,JWL,JWR,WL,JPER)
      ! ST(2-7) Drl(12)-C"(22)-
      if (JB_SYS > 0) then
        IWDS = JUST(LRJ,LRI)
        IWDT = JUST(LRI,LRJ)
        VLOP1 = W1*W1_ST(7)
        VLOP0 = Zero             !D2-6
        !LIST = LIST3(LRA,LRB,LRI)
        !WL = WL-VLOP1*vint_ci(LIST)    !4.3 Vlop0=0        !!!!!
        if (LINE == 26) then    !LRI,LRA
          call COMP_LOOP_G(9,LRI,0,0,LRA,VLOP0,VLOP1,WL0,LIST0,WL1,LIST1)
        end if
        if (LINE == 28) then    !LRI,LRS,LRA
          call COMP_LOOP_G(12,LRI,0,LRS,LRA,VLOP0,VLOP1,WL0,LIST0,WL1,LIST1)
        end if
        if (LINE == 29) then    !LRI,LRS,LRA
          call COMP_LOOP_G(11,LRI,0,LRS,LRA,VLOP0,VLOP1,WL0,LIST0,WL1,LIST1)
        end if
        !call PRODAB_1(3,JPEL,IWDS,IWDT,JWL,JWR,WL2,JPER,LIST2,LIST3)
        call PRODAB_2(3,JPEL,IWDS,IWDT,JWL,JWR,WL0,JPER,LIST0)
        if (LIST1 /= 0) call PRODAB_2(3,JPEL,IWDS,IWDT,JWL,JWR,WL1,JPER,LIST1)

        !call PRODAB(3,JPEL,IWDS,IWDT,JWL,JWR,WL,JPER)
      end if
    end if

    JMLR = Mul(JML,JMR)
    if (LMIJ /= JMLR) cycle
    W1ST1 = W1_ST(1)
    W1ST2 = W1_ST(2)
    W1ST3 = W1_ST(3)
    W1ST4 = W1_ST(4)
    NI = mod(LRJ-LRI,2)
    if (NI == 0) then
      W1ST1 = -W1ST1
      W1ST2 = -W1ST2
      W1ST3 = -W1ST3
      W1ST4 = -W1ST4
    end if
    if (JML == 1) then
      ! ST(2-1) Ar(02)-Bl(32)-
      IWDS = JUST(LRI,LRI)
      IWDT = JUST(LRI,LRJ)
      VLOP1 = W1*W1ST1
      VLOP0 = Zero
      !LIST = LIST4(LRI,LRJ,LRA,LRB)
      if (LINE == 26) then    !LRI,LRJ,LRA
        call COMP_LOOP_G(5,LRI,LRJ,LRS,LRA,VLOP0,VLOP1,WL0,LIST0,WL1,LIST1)
      end if
      if (LINE == 28) then    !LRI,LRJ,LRS,LRA
        call COMP_LOOP_G(7,LRI,LRJ,LRS,LRA,VLOP0,VLOP1,WL0,LIST0,WL1,LIST1)
      end if
      if (LINE == 29) then    !LRI,LRJ,LRS,LRA
        call COMP_LOOP_G(6,LRI,LRJ,LRS,LRA,VLOP0,VLOP1,WL0,LIST0,WL1,LIST1)
      end if
      !WL = -VLOP1*vint_ci(LIST)    !1.1 VLOP0=0        !!!!!
      !call PRODAB(3,JPEL,IWDS,IWDT,JWL,JWR,WL,JPER)
      !call PRODAB_1(3,JPEL,IWDS,IWDT,JWL,JWR,WL2,JPER,LIST2,LIST3)
      call PRODAB_2(3,JPEL,IWDS,IWDT,JWL,JWR,WL0,JPER,LIST0)
      if (LIST1 /= 0) call PRODAB_2(3,JPEL,IWDS,IWDT,JWL,JWR,WL1,JPER,LIST1)

    end if

    ! ST(2-2) (22)Ar(13)-Bl(32)-
    VLOP1 = W1*W1ST2
    !WL = -VLOP1*vint_ci(LIST)    !1.1
    VLOP0 = Zero
    if (LINE == 26) then    !LRI,LRJ,LRA
      call COMP_LOOP_G(5,LRI,LRJ,LRS,LRA,VLOP0,VLOP1,WL0,LIST0,WL1,LIST1)
    end if
    if (LINE == 28) then    !LRI,LRJ,LRS,LRA
      call COMP_LOOP_G(7,LRI,LRJ,LRS,LRA,VLOP0,VLOP1,WL0,LIST0,WL1,LIST1)
    end if
    if (LINE == 29) then    !LRI,LRJ,LRS,LRA
      call COMP_LOOP_G(6,LRI,LRJ,LRS,LRA,VLOP0,VLOP1,WL0,LIST0,WL1,LIST1)
    end if
    do LRK=NORB_FRZ+1,LRI-1
      LMK = LSM_INN(LRK)
      if (Mul(LMK,LMI) /= JML) cycle
      if (Mul(LMK,LMJ) /= JMR) cycle
      IWDS = JUST(LRK,LRI)
      IWDT = JUST(LRK,LRJ)
      !call PRODAB(3,JPEL,IWDS,IWDT,JWL,JWR,WL,JPER)
      !call PRODAB_1(3,JPEL,IWDS,IWDT,JWL,JWR,WL2,JPER,LIST2,LIST3)
      call PRODAB_2(3,JPEL,IWDS,IWDT,JWL,JWR,WL0,JPER,LIST0)
      if (LIST1 /= 0) call PRODAB_2(3,JPEL,IWDS,IWDT,JWL,JWR,WL1,JPER,LIST1)
    end do
    ! ST(2-3) Ar(13)-Bl(32)-C'(22)-
    ! ST(2-4) Ar(23)-Bl(32)-C'(12)-
    VLOP1 = W1*W1ST4
    VLOP0 = Zero
    if (LINE == 26) then    !LRI,LRJ,LRA
      call COMP_LOOP_G(5,LRI,LRJ,LRS,LRA,VLOP0,VLOP1,WL0,LIST0,WL1,LIST1)
    end if
    if (LINE == 28) then    !LRI,LRJ,LRS,LRA
      call COMP_LOOP_G(7,LRI,LRJ,LRS,LRA,VLOP0,VLOP1,WL0,LIST0,WL1,LIST1)
    end if
    if (LINE == 29) then    !LRI,LRJ,LRS,LRA
      call COMP_LOOP_G(6,LRI,LRJ,LRS,LRA,VLOP0,VLOP1,WL0,LIST0,WL1,LIST1)
    end if
    !WL = -VLOP1*vint_ci(LIST)    !1.1 VLOP0=0
    if (JB_SYS > 0) then
      VLOP1 = W1*W1ST3
      VLOP0 = Zero
      if (LINE == 26) then    !LRI,LRJ,LRA
        call COMP_LOOP_G(5,LRI,LRJ,LRS,LRA,VLOP0,VLOP1,WL0_1,LIST0_1,WL1_1,LIST1_1)
      end if
      if (LINE == 28) then    !LRI,LRJ,LRS,LRA
        call COMP_LOOP_G(7,LRI,LRJ,LRS,LRA,VLOP0,VLOP1,WL0_1,LIST0_1,WL1_1,LIST1_1)
      end if
      if (LINE == 29) then    !LRI,LRJ,LRS,LRA
        call COMP_LOOP_G(6,LRI,LRJ,LRS,LRA,VLOP0,VLOP1,WL0_1,LIST0_1,WL1_1,LIST1_1)
      end if
    end if
    !WL1 = -VLOP1*vint_ci(LIST)
    do LRK=LRJ+1,NORB_DZ
      LMK = LSM_INN(LRK)
      if (Mul(LMK,LMI) /= JML) cycle
      if (Mul(LMK,LMJ) /= JMR) cycle
      IWDS = JUST(LRI,LRK)
      IWDT = JUST(LRJ,LRK)
      !call PRODAB(3,JPEL,IWDS,IWDT,JWL,JWR,WL,JPER)
      !call PRODAB_1(3,JPEL,IWDS,IWDT,JWL,JWR,WL2,JPER,LIST2,LIST3)
      call PRODAB_2(3,JPEL,IWDS,IWDT,JWL,JWR,WL0,JPER,LIST0)
      if (LIST1 /= 0) call PRODAB_2(3,JPEL,IWDS,IWDT,JWL,JWR,WL1,JPER,LIST1)

      if (JB_SYS > 0) then
        IWDS = JUST(LRK,LRI)
        !call PRODAB(3,JPEL,IWDS,IWDT,JWL,JWR,WL1,JPER)
        !call PRODAB_1(3,JPEL,IWDS,IWDT,JWL,JWR,WL2_1,JPER,LIST2_1,LIST3_1)
        call PRODAB_2(3,JPEL,IWDS,IWDT,JWL,JWR,WL0_1,JPER,LIST0_1)
        if (LIST1_1 /= 0) call PRODAB_2(3,JPEL,IWDS,IWDT,JWL,JWR,WL1_1,JPER,LIST1_1)

      end if
    end do
    ! ST(2-4) Ar(23)-C'(12)-Bl(32)-
    do LRK=LRI+1,LRJ-1
      LMK = LSM_INN(LRK)
      if (Mul(LMK,LMI) /= JML) cycle
      if (Mul(LMK,LMJ) /= JMR) cycle
      IWDS = JUST(LRI,LRK)
      IWDT = JUST(LRK,LRJ)
      !call PRODAB(3,JPEL,IWDS,IWDT,JWL,JWR,-WL,JPER)
      !call PRODAB_1(3,JPEL,IWDS,IWDT,JWL,JWR,-WL2,JPER,LIST2,LIST3
      call PRODAB_2(3,JPEL,IWDS,IWDT,JWL,JWR,-WL0,JPER,LIST0)
      if (LIST1 /= 0) call PRODAB_2(3,JPEL,IWDS,IWDT,JWL,JWR,-WL1,JPER,LIST1)

      ! ST(2-3) Ar(13)-C'(22)-Bl(32)-
      if (JB_SYS > 0) then
        IWDS = JUST(LRK,LRI)
        !call PRODAB(3,JPEL,IWDS,IWDT,JWL,JWR,-WL1,JPER)
        !call PRODAB_1(3,JPEL,IWDS,IWDT,JWL,JWR,-WL2_1,JPER,LIST2_1,LIST3_1)
        call PRODAB_2(3,JPEL,IWDS,IWDT,JWL,JWR,-WL0_1,JPER,LIST0_1)
        if (LIST1_1 /= 0) call PRODAB_2(3,JPEL,IWDS,IWDT,JWL,JWR,-WL1_1,JPER,LIST1_1)

      end if
    end do
  end do
end do

return

end subroutine ST_HEAD_DBL_TAIL_ACT_G

subroutine TS_HEAD_DBL_TAIL_ACT_G(LRA)

use gugaci_global, only: jb_sys, jml, jmr, jpel, jper, just, jwl, jwr, line, lrs, lsm_inn, norb_dz, norb_frz, w1, w1_ts
use Symmetry_Info, only: Mul
use Constants, only: Zero
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: lra
integer(kind=iwp) :: iwds, iwdt, list0, list1, lmi, lmij, lmj, lmk, lri, lrj, lrk, ni
real(kind=wp) :: vlop0, vlop1, w1ts1, w1ts2, w1ts3, w1ts4, wl0, wl1

! TS(3-1) Ar(23)-Bl(20)-
! TS(3-2) (22)Ar(23)-Bl(31)-
! TS(3-2) Ar(23)-C'(22)-Bl(31)-
! TS(3-3) Ar(23)-Bl(31)-C"(22)-
! TS(3-4) Ar(23)-Bl(32)-C"(21)-
do LRI=NORB_FRZ+1,NORB_DZ
  LMI = LSM_INN(LRI)
  do LRJ=LRI+1,NORB_DZ
    LMJ = LSM_INN(LRJ)
    LMIJ = Mul(LMI,LMJ)
    W1TS1 = W1_TS(1)
    W1TS2 = W1_TS(2)
    W1TS3 = W1_TS(3)
    W1TS4 = W1_TS(4)
    NI = mod(LRJ-LRI,2)
    if (NI == 0) then
      W1TS1 = -W1TS1
      W1TS2 = -W1TS2
      W1TS3 = -W1TS3
      W1TS4 = -W1TS4
    end if
    !LIST = LIST3(LRI,LRJ,LRB)
    !-------------------------------------------------------------------
    ! TS(3-1) Ar(23)-Bl(20)-
    if ((JMR == 1) .and. (LMIJ == JML)) then
      IWDT = JUST(LRI,LRJ)
      IWDS = JUST(LRJ,LRJ)
      VLOP1 = W1*W1TS1             !D3-1
      VLOP0 = Zero
      if (LINE == 26) then    !LRI,LRJ,LRA
        call COMP_LOOP_G(5,LRI,LRJ,LRS,LRA,VLOP0,VLOP1,WL0,LIST0,WL1,LIST1)
      end if
      if (LINE == 28) then    !LRI,LRJ,LRS,LRA
        call COMP_LOOP_G(7,LRI,LRJ,LRS,LRA,VLOP0,VLOP1,WL0,LIST0,WL1,LIST1)
      end if
      if (LINE == 29) then    !LRI,LRJ,LRS,LRA
        call COMP_LOOP_G(6,LRI,LRJ,LRS,LRA,VLOP0,VLOP1,WL0,LIST0,WL1,LIST1)
      end if
      !WL=  -VLOP1*(VINT_CI(LIST))    !2.2 vlop0=0
      !call PRODAB(3,JPEL,IWDT,IWDS,JWL,JWR,WL,JPER)
      !call PRODAB_1(3,JPEL,IWDT,IWDS,JWL,JWR,WL2,JPER,LIST2,LIST3)
      call PRODAB_2(3,JPEL,IWDT,IWDS,JWL,JWR,WL0,JPER,LIST0)
      if (LIST1 /= 0) call PRODAB_2(3,JPEL,IWDT,IWDS,JWL,JWR,WL1,JPER,LIST1)
    end if
    !-------------------------------------------------------------------
    VLOP1 = W1*W1TS2             !D3-2
    VLOP0 = Zero
    if (LINE == 26) then    !LRI,LRJ,LRA
      call COMP_LOOP_G(5,LRI,LRJ,LRS,LRA,VLOP0,VLOP1,WL0,LIST0,WL1,LIST1)
    end if
    if (LINE == 28) then    !LRI,LRJ,LRS,LRA
      call COMP_LOOP_G(7,LRI,LRJ,LRS,LRA,VLOP0,VLOP1,WL0,LIST0,WL1,LIST1)
    end if
    if (LINE == 29) then    !LRI,LRJ,LRS,LRA
      call COMP_LOOP_G(6,LRI,LRJ,LRS,LRA,VLOP0,VLOP1,WL0,LIST0,WL1,LIST1)
    end if
    !WL = -VLOP1*(VINT_CI(LIST))    !2.2 vlop0=0
    ! TS(3-2) (22)Ar(23)-Bl(31)-
    do LRK=NORB_FRZ+1,LRI-1
      LMK = LSM_INN(LRK)
      if ((Mul(LMK,LMI) /= JML) .or. (Mul(LMK,LMJ) /= JMR)) cycle
      IWDT = JUST(LRK,LRI)
      IWDS = JUST(LRK,LRJ)
      !call PRODAB(3,JPEL,IWDT,IWDS,JWL,JWR,WL,JPER)
      !call PRODAB_1(3,JPEL,IWDT,IWDS,JWL,JWR,WL2,JPER,LIST2,LIST3)
      call PRODAB_2(3,JPEL,IWDT,IWDS,JWL,JWR,WL0,JPER,LIST0)
      if (LIST1 /= 0) call PRODAB_2(3,JPEL,IWDT,IWDS,JWL,JWR,WL1,JPER,LIST1)
    end do

    ! TS(3-2) Ar(23)-C'(22)-Bl(31)-
    do LRK=LRI+1,LRJ-1
      LMK = LSM_INN(LRK)
      if ((Mul(LMK,LMI) /= JML) .or. (Mul(LMK,LMJ) /= JMR)) cycle
      IWDT = JUST(LRI,LRK)
      IWDS = JUST(LRK,LRJ)
      !call PRODAB(3,JPEL,IWDT,IWDS,JWL,JWR,-WL,JPER)
      !call PRODAB_1(3,JPEL,IWDT,IWDS,JWL,JWR,-WL2,JPER,LIST2,LIST3
      call PRODAB_2(3,JPEL,IWDT,IWDS,JWL,JWR,-WL0,JPER,LIST0)
      if (LIST1 /= 0) call PRODAB_2(3,JPEL,IWDT,IWDS,JWL,JWR,-WL1,JPER,LIST1)
    end do

    !-------------------------------------------------------------------
    ! TS(3-4) Ar(23)-Bl(32)-C"(21)-
    VLOP1 = W1*W1TS4             !D3-4
    VLOP0 = Zero
    if (LINE == 26) then    !LRI,LRJ,LRA
      call COMP_LOOP_G(5,LRI,LRJ,LRS,LRA,VLOP0,VLOP1,WL0,LIST0,WL1,LIST1)
    end if
    if (LINE == 28) then    !LRI,LRJ,LRS,LRA
      call COMP_LOOP_G(7,LRI,LRJ,LRS,LRA,VLOP0,VLOP1,WL0,LIST0,WL1,LIST1)
    end if
    if (LINE == 29) then    !LRI,LRJ,LRS,LRA
      call COMP_LOOP_G(6,LRI,LRJ,LRS,LRA,VLOP0,VLOP1,WL0,LIST0,WL1,LIST1)
    end if
    do LRK=LRJ+1,NORB_DZ
      LMK = LSM_INN(LRK)
      if ((Mul(LMK,LMI) /= JML) .or. (Mul(LMK,LMJ) /= JMR)) cycle
      IWDT = JUST(LRI,LRK)
      IWDS = JUST(LRJ,LRK)
      !call PRODAB(3,JPEL,IWDT,IWDS,JWL,JWR,WL,JPER)
      !call PRODAB_1(3,JPEL,IWDT,IWDS,JWL,JWR,WL2,JPER,LIST2,LIST3)
      call PRODAB_2(3,JPEL,IWDT,IWDS,JWL,JWR,WL0,JPER,LIST0)
      if (LIST1 /= 0) call PRODAB_2(3,JPEL,IWDT,IWDS,JWL,JWR,WL1,JPER,LIST1)

    end do

    !-------------------------------------------------------------------
    ! TS(3-3) Ar(23)-Bl(31)-C"(22)-
    if (JB_SYS > 0) then
      VLOP1 = W1*W1TS3             !D3-4
      VLOP0 = Zero
      if (LINE == 26) then    !LRI,LRJ,LRA
        call COMP_LOOP_G(5,LRI,LRJ,LRS,LRA,VLOP0,VLOP1,WL0,LIST0,WL1,LIST1)
      end if
      if (LINE == 28) then    !LRI,LRJ,LRS,LRA
        call COMP_LOOP_G(7,LRI,LRJ,LRS,LRA,VLOP0,VLOP1,WL0,LIST0,WL1,LIST1)
      end if
      if (LINE == 29) then    !LRI,LRJ,LRS,LRA
        call COMP_LOOP_G(6,LRI,LRJ,LRS,LRA,VLOP0,VLOP1,WL0,LIST0,WL1,LIST1)
      end if
      !WL = -VLOP1*VINT_CI(LIST)    !2.2 vlop0=0
      do LRK=LRJ+1,NORB_DZ
        LMK = LSM_INN(LRK)
        if ((Mul(LMK,LMI) /= JML) .or. (Mul(LMK,LMJ) /= JMR)) cycle
        IWDT = JUST(LRI,LRK)
        IWDS = JUST(LRK,LRJ)
        !call PRODAB(3,JPEL,IWDT,IWDS,JWL,JWR,WL,JPER)
        !call PRODAB_1(3,JPEL,IWDT,IWDS,JWL,JWR,WL2,JPER,LIST2,LIST3)
        call PRODAB_2(3,JPEL,IWDT,IWDS,JWL,JWR,WL0,JPER,LIST0)
        if (LIST1 /= 0) call PRODAB_2(3,JPEL,IWDT,IWDS,JWL,JWR,WL1,JPER,LIST1)

      end do
    end if
  end do
end do

return

end subroutine TS_HEAD_DBL_TAIL_ACT_G

subroutine STT_HEAD_DBL_TAIL_ACT_G(LRA)

use gugaci_global, only: jml, jmr, jpel, jper, just, jwl, jwr, line, lrs, lsm_inn, norb_dz, norb_frz, w1, w1_st1
use Symmetry_Info, only: Mul
use Constants, only: Zero
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: lra
integer(kind=iwp) :: iwdl, iwdr, list0, list1, lmi, lmij, lmj, lmk, lri, lrj, lrk, ni
real(kind=wp) :: vlop0, vlop1, w1st1, w1st2, w1st3, w1st4, wl0, wl1

! ST1(4-1) Ar(01)-Bl(31)-
! ST1(4-2) (11)Ar(23)-Bl(31)-
! ST1(4-3) Ar(13)-C'(21)-Bl(31)-
! ST1(4-3) Ar(13)-Bl(31)-C"(21)-
! ST1(4-4) Ar(23)-C'(11)-Bl(31)-
! ST1(4-4) Ar(23)-Bl(31)-C"(11)-
do LRI=NORB_FRZ+1,NORB_DZ-1
  LMI = LSM_INN(LRI)
  do LRJ=LRI+1,NORB_DZ
    LMJ = LSM_INN(LRJ)
    LMIJ = Mul(LMI,LMJ)
    W1ST1 = W1_ST1(1)
    W1ST2 = W1_ST1(2)
    W1ST3 = W1_ST1(3)
    W1ST4 = W1_ST1(4)
    NI = mod(LRJ-LRI,2)
    if (NI == 0) then
      W1ST1 = -W1ST1
      W1ST2 = -W1ST2
      W1ST3 = -W1ST3
      W1ST4 = -W1ST4
    end if
    !      LIST=LIST3(LRI,LRJ,LRB)
    if ((JML == 1) .and. (LMIJ == JMR)) then
      ! ST1(4-1) Ar(01)-Bl(31)-
      IWDL = JUST(LRI,LRI)
      IWDR = JUST(LRI,LRJ)
      VLOP1 = W1*W1ST1
      VLOP0 = Zero
      if (LINE == 26) then    !LRI,LRJ,LRA
        call COMP_LOOP_G(5,LRI,LRJ,LRS,LRA,VLOP0,VLOP1,WL0,LIST0,WL1,LIST1)
      end if
      if (LINE == 28) then    !LRI,LRJ,LRS,LRA
        call COMP_LOOP_G(7,LRI,LRJ,LRS,LRA,VLOP0,VLOP1,WL0,LIST0,WL1,LIST1)
      end if
      if (LINE == 29) then    !LRI,LRJ,LRS,LRA
        call COMP_LOOP_G(6,LRI,LRJ,LRS,LRA,VLOP0,VLOP1,WL0,LIST0,WL1,LIST1)
      end if
      !WL = -VLOP1*(VINT_CI(LIST))    !2.2 vlop0=0
      !call PRODAB(3,JPEL,IWDL,IWDR,JWL,JWR,WL,JPER)
      !call PRODAB_1(3,JPEL,IWDL,IWDR,JWL,JWR,WL2,JPER,LIST2,LIST3)
      call PRODAB_2(3,JPEL,IWDL,IWDR,JWL,JWR,WL0,JPER,LIST0)
      if (LIST1 /= 0) call PRODAB_2(3,JPEL,IWDL,IWDR,JWL,JWR,WL1,JPER,LIST1)
    end if
    ! ST1(4-2) (11)Ar(23)-Bl(31)-
    VLOP1 = W1*W1ST2
    VLOP0 = Zero
    if (LINE == 26) then    !LRI,LRJ,LRA
      call COMP_LOOP_G(5,LRI,LRJ,LRS,LRA,VLOP0,VLOP1,WL0,LIST0,WL1,LIST1)
    end if
    if (LINE == 28) then    !LRI,LRJ,LRS,LRA
      call COMP_LOOP_G(7,LRI,LRJ,LRS,LRA,VLOP0,VLOP1,WL0,LIST0,WL1,LIST1)
    end if
    if (LINE == 29) then    !LRI,LRJ,LRS,LRA
      call COMP_LOOP_G(6,LRI,LRJ,LRS,LRA,VLOP0,VLOP1,WL0,LIST0,WL1,LIST1)
    end if
    do LRK=NORB_FRZ+1,LRI-1
      LMK = LSM_INN(LRK)
      if ((Mul(LMK,LMI) /= JML) .or. (Mul(LMK,LMJ) /= JMR)) cycle
      !WL = -VLOP1*VINT_CI(LIST)
      IWDL = JUST(LRI,LRK)
      IWDR = JUST(LRK,LRJ)
      !call PRODAB(3,JPEL,IWDL,IWDR,JWL,JWR,WL,JPER)
      !call PRODAB_1(3,JPEL,IWDL,IWDR,JWL,JWR,WL2,JPER,LIST2,LIST3)
      call PRODAB_2(3,JPEL,IWDL,IWDR,JWL,JWR,WL0,JPER,LIST0)
      if (LIST1 /= 0) call PRODAB_2(3,JPEL,IWDL,IWDR,JWL,JWR,WL1,JPER,LIST1)
    end do
    VLOP1 = W1*W1ST4
    VLOP0 = Zero
    if (LINE == 26) then    !LRI,LRJ,LRA
      call COMP_LOOP_G(5,LRI,LRJ,LRS,LRA,VLOP0,VLOP1,WL0,LIST0,WL1,LIST1)
    end if
    if (LINE == 28) then    !LRI,LRJ,LRS,LRA
      call COMP_LOOP_G(7,LRI,LRJ,LRS,LRA,VLOP0,VLOP1,WL0,LIST0,WL1,LIST1)
    end if
    if (LINE == 29) then    !LRI,LRJ,LRS,LRA
      call COMP_LOOP_G(6,LRI,LRJ,LRS,LRA,VLOP0,VLOP1,WL0,LIST0,WL1,LIST1)
    end if
    ! ST1(4-4) Ar(23)-C'(11)-Bl(31)-
    do LRK=LRI+1,LRJ-1
      LMK = LSM_INN(LRK)
      if ((Mul(LMI,LMK) /= JML) .or. (Mul(LMK,LMJ) /= JMR)) cycle
      IWDL = JUST(LRI,LRK)
      IWDR = JUST(LRK,LRJ)
      !call PRODAB(3,JPEL,IWDL,IWDR,JWL,JWR,-WL,JPER)
      !call PRODAB_1(3,JPEL,IWDL,IWDR,JWL,JWR,-WL2,JPER,LIST2,LIST3
      call PRODAB_2(3,JPEL,IWDL,IWDR,JWL,JWR,-WL0,JPER,LIST0)
      if (LIST1 /= 0) call PRODAB_2(3,JPEL,IWDL,IWDR,JWL,JWR,-WL1,JPER,LIST1)
    end do
    ! ST1(4-4) Ar(23)-Bl(31)-C"(11)-
    do LRK=LRJ+1,NORB_DZ
      LMK = LSM_INN(LRK)
      if ((Mul(LMI,LMK) /= JML) .or. (Mul(LMJ,LMK) /= JMR)) cycle
      IWDL = JUST(LRI,LRK)
      IWDR = JUST(LRJ,LRK)
      !call PRODAB(3,JPEL,IWDL,IWDR,JWL,JWR,WL,JPER)
      !call PRODAB_1(3,JPEL,IWDL,IWDR,JWL,JWR,WL2,JPER,LIST2,LIST3)
      call PRODAB_2(3,JPEL,IWDL,IWDR,JWL,JWR,WL0,JPER,LIST0)
      if (LIST1 /= 0) call PRODAB_2(3,JPEL,IWDL,IWDR,JWL,JWR,WL1,JPER,LIST1)
    end do

    VLOP1 = W1*W1ST3
    VLOP0 = Zero
    if (LINE == 26) then    !LRI,LRJ,LRA
      call COMP_LOOP_G(5,LRI,LRJ,LRS,LRA,VLOP0,VLOP1,WL0,LIST0,WL1,LIST1)
    end if
    if (LINE == 28) then    !LRI,LRJ,LRS,LRA
      call COMP_LOOP_G(7,LRI,LRJ,LRS,LRA,VLOP0,VLOP1,WL0,LIST0,WL1,LIST1)
    end if
    if (LINE == 29) then    !LRI,LRJ,LRS,LRA
      call COMP_LOOP_G(6,LRI,LRJ,LRS,LRA,VLOP0,VLOP1,WL0,LIST0,WL1,LIST1)
    end if
    ! ST1(4-3) Ar(13)-C'(21)-Bl(31)-
    do LRK=LRI+1,LRJ-1
      LMK = LSM_INN(LRK)
      if ((Mul(LMI,LMK) /= JML) .or. (Mul(LMK,LMJ) /= JMR)) cycle
      IWDL = JUST(LRK,LRI)
      IWDR = JUST(LRK,LRJ)
      !call PRODAB(3,JPEL,IWDL,IWDR,JWL,JWR,-WL,JPER)
      !call PRODAB_1(3,JPEL,IWDL,IWDR,JWL,JWR,-WL2,JPER,LIST2,LIST3
      call PRODAB_2(3,JPEL,IWDL,IWDR,JWL,JWR,-WL0,JPER,LIST0)
      if (LIST1 /= 0) call PRODAB_2(3,JPEL,IWDL,IWDR,JWL,JWR,-WL1,JPER,LIST1)
    end do
    ! ST1(4-3) Ar(13)-Bl(31)-C"(21)-
    do LRK=LRJ+1,NORB_DZ
      LMK = LSM_INN(LRK)
      if ((Mul(LMI,LMK) /= JML) .or. (Mul(LMJ,LMK) /= JMR)) cycle
      IWDL = JUST(LRK,LRI)
      IWDR = JUST(LRJ,LRK)
      !call PRODAB(3,JPEL,IWDL,IWDR,JWL,JWR,WL,JPER)
      !call PRODAB_1(3,JPEL,IWDL,IWDR,JWL,JWR,WL2,JPER,LIST2,LIST3)
      call PRODAB_2(3,JPEL,IWDL,IWDR,JWL,JWR,WL0,JPER,LIST0)
      if (LIST1 /= 0) call PRODAB_2(3,JPEL,IWDL,IWDR,JWL,JWR,WL1,JPER,LIST1)
    end do
  end do
end do

return

end subroutine STT_HEAD_DBL_TAIL_ACT_G

subroutine TTS_HEAD_DBL_TAIL_ACT_G(LRA)

use gugaci_global, only: jml, jmr, jpel, jper, just, jwl, jwr, line, lrs, lsm_inn, norb_dz, norb_frz, w1, w1_t1s
use Symmetry_Info, only: Mul
use Constants, only: Zero
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: lra
integer(kind=iwp) :: iwdl, iwdr, list0, list1, lmi, lmij, lmj, lmk, lri, lrj, lrk, ni
real(kind=wp) :: vlop0, vlop1, w1t1s1, w1t1s2, w1t1s3, w1t1s4, w1t1s5, w1t1s6, w1t1s7, wl0, wl1

! T1S(5-1)   Ar(13)-Bl(10)-
! T1S(5-2)   Ar(13)-Bl(32)-
! T1S(5-2)   Ar(13)-C'(11)-Bl(32)-
! T1S(5-3)   Ar(13)-Bl(31)-C"(12)-
! T1S(5-4)   Ar(13)-Bl(32)-C"(11)-
! T1S(5-5)   Drl(12)-
! T1S(5-6)   Drl(12)-C"(12)-
! T1S(5-7)   Drl(12)-C"(11)-
do LRI=NORB_FRZ+1,NORB_DZ
  LMI = LSM_INN(LRI)
  do LRJ=LRI+1,NORB_DZ
    LMJ = LSM_INN(LRJ)
    LMIJ = Mul(LMI,LMJ)
    W1T1S1 = W1_T1S(1)
    W1T1S2 = W1_T1S(2)
    W1T1S3 = W1_T1S(3)
    W1T1S4 = W1_T1S(4)
    W1T1S5 = W1_T1S(5)
    W1T1S6 = W1_T1S(6)
    W1T1S7 = W1_T1S(7)
    NI = mod(LRJ-LRI,2)
    if (NI == 0) then
      W1T1S1 = -W1T1S1
      W1T1S2 = -W1T1S2
      W1T1S3 = -W1T1S3
      W1T1S4 = -W1T1S4
      !W1T1S5 = -W1T1S5
      !W1T1S6 = -W1T1S6
      !W1T1S7 = -W1T1S7
    end if
    !LIST = LIST3(LRI,LRJ,LRB)
    if ((JML == LMIJ) .and. (JMR == 1)) then
      ! T1S(5-1)   Ar(13)-Bl(10)-
      VLOP1 = W1*W1T1S1
      VLOP0 = Zero
      if (LINE == 26) then    !LRI,LRJ,LRA
        call COMP_LOOP_G(5,LRI,LRJ,LRS,LRA,VLOP0,VLOP1,WL0,LIST0,WL1,LIST1)
      end if
      if (LINE == 28) then    !LRI,LRJ,LRS,LRA
        call COMP_LOOP_G(7,LRI,LRJ,LRS,LRA,VLOP0,VLOP1,WL0,LIST0,WL1,LIST1)
      end if
      if (LINE == 29) then    !LRI,LRJ,LRS,LRA
        call COMP_LOOP_G(6,LRI,LRJ,LRS,LRA,VLOP0,VLOP1,WL0,LIST0,WL1,LIST1)
      end if
      !WL = -VLOP1*VINT_CI(LIST)
      IWDL = JUST(LRI,LRJ)
      IWDR = JUST(LRJ,LRJ)
      !call PRODAB(3,JPEL,IWDL,IWDR,JWL,JWR,WL,JPER)
      !call PRODAB_1(3,JPEL,IWDL,IWDR,JWL,JWR,WL2,JPER,LIST2,LIST3)
      call PRODAB_2(3,JPEL,IWDL,IWDR,JWL,JWR,WL0,JPER,LIST0)
      if (LIST1 /= 0) call PRODAB_2(3,JPEL,IWDL,IWDR,JWL,JWR,WL1,JPER,LIST1)
    end if
    VLOP1 = W1*W1T1S2
    VLOP0 = Zero
    if (LINE == 26) then    !LRI,LRJ,LRA
      call COMP_LOOP_G(5,LRI,LRJ,LRS,LRA,VLOP0,VLOP1,WL0,LIST0,WL1,LIST1)
    end if
    if (LINE == 28) then    !LRI,LRJ,LRS,LRA
      call COMP_LOOP_G(7,LRI,LRJ,LRS,LRA,VLOP0,VLOP1,WL0,LIST0,WL1,LIST1)
    end if
    if (LINE == 29) then    !LRI,LRJ,LRS,LRA
      call COMP_LOOP_G(6,LRI,LRJ,LRS,LRA,VLOP0,VLOP1,WL0,LIST0,WL1,LIST1)
    end if
    !WL = -VLOP1*VINT_CI(LIST)
    do LRK=NORB_FRZ+1,LRI-1
      LMK = LSM_INN(LRK)
      ! T1S(5-2)   (11)Ar(13)-Bl(32)-
      if ((JML /= Mul(LMK,LMI)) .or. (JMR /= Mul(LMK,LMJ))) cycle
      !VLOP1 = W1*W1T1S2
      !WL = -VLOP1*VINT_CI(LIST)
      IWDL = JUST(LRK,LRI)
      IWDR = JUST(LRJ,LRK)
      !call PRODAB(3,JPEL,IWDL,IWDR,JWL,JWR,WL,JPER)
      !call PRODAB_1(3,JPEL,IWDL,IWDR,JWL,JWR,WL2,JPER,LIST2,LIST3)
      call PRODAB_2(3,JPEL,IWDL,IWDR,JWL,JWR,WL0,JPER,LIST0)
      if (LIST1 /= 0) call PRODAB_2(3,JPEL,IWDL,IWDR,JWL,JWR,WL1,JPER,LIST1)

    end do
    do LRK=LRI+1,LRJ-1
      LMK = LSM_INN(LRK)
      ! T1S(5-2)   Ar(13)-C'(11)-Bl(32)-
      if ((JML /= Mul(LMI,LMK)) .or. (JMR /= Mul(LMK,LMJ))) cycle
      !VLOP1 = W1*W1T1S2
      !WL = -VLOP1*VINT_CI(LIST)
      IWDL = JUST(LRI,LRK)
      IWDR = JUST(LRJ,LRK)
      !call PRODAB(3,JPEL,IWDL,IWDR,JWL,JWR,-WL,JPER)
      !call PRODAB_1(3,JPEL,IWDL,IWDR,JWL,JWR,-WL2,JPER,LIST2,LIST3
      call PRODAB_2(3,JPEL,IWDL,IWDR,JWL,JWR,-WL0,JPER,LIST0)
      if (LIST1 /= 0) call PRODAB_2(3,JPEL,IWDL,IWDR,JWL,JWR,-WL1,JPER,LIST1)

    end do
    VLOP1 = W1*W1T1S3
    VLOP0 = Zero
    if (LINE == 26) then    !LRI,LRJ,LRA
      call COMP_LOOP_G(5,LRI,LRJ,LRS,LRA,VLOP0,VLOP1,WL0,LIST0,WL1,LIST1)
    end if
    if (LINE == 28) then    !LRI,LRJ,LRS,LRA
      call COMP_LOOP_G(7,LRI,LRJ,LRS,LRA,VLOP0,VLOP1,WL0,LIST0,WL1,LIST1)
    end if
    if (LINE == 29) then    !LRI,LRJ,LRS,LRA
      call COMP_LOOP_G(6,LRI,LRJ,LRS,LRA,VLOP0,VLOP1,WL0,LIST0,WL1,LIST1)
    end if
    do LRK=LRJ+1,NORB_DZ
      LMK = LSM_INN(LRK)
      ! T1S(5-3)   Ar(13)-Bl(31)-C"(12)-
      if ((JML /= Mul(LMI,LMK)) .or. (JMR /= Mul(LMJ,LMK))) cycle
      !VLOP1 = W1*W1T1S3
      !WL = -VLOP1*VINT_CI(LIST)
      IWDL = JUST(LRI,LRK)
      IWDR = JUST(LRK,LRJ)
      !call PRODAB(3,JPEL,IWDL,IWDR,JWL,JWR,WL,JPER)
      !call PRODAB_1(3,JPEL,IWDL,IWDR,JWL,JWR,WL2,JPER,LIST2,LIST3)
      call PRODAB_2(3,JPEL,IWDL,IWDR,JWL,JWR,WL0,JPER,LIST0)
      if (LIST1 /= 0) call PRODAB_2(3,JPEL,IWDL,IWDR,JWL,JWR,WL1,JPER,LIST1)

    end do
    VLOP1 = W1*W1T1S4
    VLOP0 = Zero
    if (LINE == 26) then    !LRI,LRJ,LRA
      call COMP_LOOP_G(5,LRI,LRJ,LRS,LRA,VLOP0,VLOP1,WL0,LIST0,WL1,LIST1)
    end if
    if (LINE == 28) then    !LRI,LRK,LRS,LRA
      call COMP_LOOP_G(7,LRI,LRJ,LRS,LRA,VLOP0,VLOP1,WL0,LIST0,WL1,LIST1)
    end if
    if (LINE == 29) then    !LRI,LRK,LRS,LRA
      call COMP_LOOP_G(6,LRI,LRJ,LRS,LRA,VLOP0,VLOP1,WL0,LIST0,WL1,LIST1)
    end if
    !WL = -VLOP1*VINT_CI(LIST)
    do LRK=LRJ+1,NORB_DZ
      LMK = LSM_INN(LRK)
      ! T1S(5-4)   Ar(13)-Bl(32)-C"(11)-
      if ((JML /= Mul(LMI,LMK)) .or. (JMR /= Mul(LMJ,LMK))) cycle
      IWDL = JUST(LRI,LRK)
      IWDR = JUST(LRJ,LRK)
      !call PRODAB(3,JPEL,IWDL,IWDR,JWL,JWR,WL,JPER)
      !call PRODAB_1(3,JPEL,IWDL,IWDR,JWL,JWR,WL2,JPER,LIST2,LIST3)
      call PRODAB_2(3,JPEL,IWDL,IWDR,JWL,JWR,WL0,JPER,LIST0)
      if (LIST1 /= 0) call PRODAB_2(3,JPEL,IWDL,IWDR,JWL,JWR,WL1,JPER,LIST1)
    end do
    if ((LMIJ == JML) .and. (LMIJ == JMR)) then
      ! T1S(5-5)   (11)Drl(12)-
      VLOP1 = W1*W1T1S5
      VLOP0 = Zero
      if (LINE == 26) then    !LRJ,LRA
        call COMP_LOOP_G(9,LRJ,0,0,LRA,VLOP0,VLOP1,WL0,LIST0,WL1,LIST1)
      end if
      if (LINE == 28) then    !LRJ,LRS,LRA
        call COMP_LOOP_G(12,LRJ,0,LRS,LRA,VLOP0,VLOP1,WL0,LIST0,WL1,LIST1)
      end if
      if (LINE == 29) then    !LRJ,LRS,LRA
        call COMP_LOOP_G(11,LRJ,0,LRS,LRA,VLOP0,VLOP1,WL0,LIST0,WL1,LIST1)
      end if
      !WL = -VLOP1*VOINT(LRB,LRJ)
      IWDL = JUST(LRI,LRJ)
      IWDR = JUST(LRJ,LRI)
      !call PRODAB(3,JPEL,IWDL,IWDR,JWL,JWR,WL,JPER)
      !call PRODAB_1(3,JPEL,IWDL,IWDR,JWL,JWR,WL2,JPER,LIST2,LIST3)
      call PRODAB_2(3,JPEL,IWDL,IWDR,JWL,JWR,WL0,JPER,LIST0)
      if (LIST1 /= 0) call PRODAB_2(3,JPEL,IWDL,IWDR,JWL,JWR,WL1,JPER,LIST1)

      ! T1S(5-6)   Drl(11)-C"(12)-
      VLOP1 = W1*W1T1S6
      VLOP0 = Zero
      if (LINE == 26) then    !LRI,LRA
        call COMP_LOOP_G(9,LRI,0,0,LRA,VLOP0,VLOP1,WL0,LIST0,WL1,LIST1)
      end if
      if (LINE == 28) then    !LRI,LRS,LRA
        call COMP_LOOP_G(12,LRI,0,LRS,LRA,VLOP0,VLOP1,WL0,LIST0,WL1,LIST1)
      end if
      if (LINE == 29) then    !LRI,LRS,LRA
        call COMP_LOOP_G(11,LRI,0,LRS,LRA,VLOP0,VLOP1,WL0,LIST0,WL1,LIST1)
      end if
      !WL = -VLOP1*VOINT(LRB,LRI)
      IWDL = JUST(LRI,LRJ)
      IWDR = JUST(LRJ,LRI)
      !call PRODAB(3,JPEL,IWDL,IWDR,JWL,JWR,WL,JPER)
      !call PRODAB_1(3,JPEL,IWDL,IWDR,JWL,JWR,WL2,JPER,LIST2,LIST3)
      call PRODAB_2(3,JPEL,IWDL,IWDR,JWL,JWR,WL0,JPER,LIST0)
      if (LIST1 /= 0) call PRODAB_2(3,JPEL,IWDL,IWDR,JWL,JWR,WL1,JPER,LIST1)

      ! T1S(5-7)   Drl(12)-C"(11)-
      VLOP1 = W1*W1T1S7
      VLOP0 = Zero
      if (LINE == 26) then    !LRI,LRA
        call COMP_LOOP_G(9,LRI,0,0,LRA,VLOP0,VLOP1,WL0,LIST0,WL1,LIST1)
      end if
      if (LINE == 28) then    !LRI,LRS,LRA
        call COMP_LOOP_G(12,LRI,0,LRS,LRA,VLOP0,VLOP1,WL0,LIST0,WL1,LIST1)
      end if
      if (LINE == 29) then    !LRI,LRS,LRA
        call COMP_LOOP_G(11,LRI,0,LRS,LRA,VLOP0,VLOP1,WL0,LIST0,WL1,LIST1)
      end if
      !WL = -VLOP1*VOINT(LRB,LRI)
      IWDL = JUST(LRI,LRJ)
      IWDR = IWDL
      !call PRODAB(3,JPEL,IWDL,IWDR,JWL,JWR,WL,JPER)
      !call PRODAB_1(3,JPEL,IWDL,IWDR,JWL,JWR,WL2,JPER,LIST2,LIST3)
      call PRODAB_2(3,JPEL,IWDL,IWDR,JWL,JWR,WL0,JPER,LIST0)
      if (LIST1 /= 0) call PRODAB_2(3,JPEL,IWDL,IWDR,JWL,JWR,WL1,JPER,LIST1)

    end if
  end do
end do

return

end subroutine TTS_HEAD_DBL_TAIL_ACT_G

subroutine TT_HEAD_DBL_TAIL_ACT_G(LRA)

use gugaci_global, only: jml, jmr, jpad, jpadl, jpel, jper, just, jwl, jwr, line, lrs, lsm_inn, norb_dz, norb_frz, w0, w0_tt, w1, &
                         w1_tt
use Symmetry_Info, only: Mul
use Constants, only: Zero
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: lra
integer(kind=iwp) :: iwdl, iwdr, jmlr, list0, list1, lmi, lmij, lmj, lmk, lmki, lmkj, lri, lrj, lrk, ni
real(kind=wp) :: vlop0, vlop1, w0tt1, w0tt2, w0tt3, w1tt1, w1tt2, wl0, wl1

! TT(11-1) (22)Ar(23)-Bl(32)-
! TT(11-1) Ar(23)-C'(22)-Bl(32)-
! TT(11-1) Ar(23)-Bl(32)-C"(22)-
! TT(11-2) (22)Drl(22)-
! TT(11-2) Drl(22)-C"(22)-
! TT(11-3) (22)Drl(33)-
! TT(11-3) Drl(33)-C"(22)-
! TT(11-3) Drl(33)-C"(22)-C"(22)-
do LRI=NORB_FRZ+1,NORB_DZ-1
  LMI = LSM_INN(LRI)
  do LRJ=LRI+1,NORB_DZ
    LMJ = LSM_INN(LRJ)
    LMIJ = Mul(LMI,LMJ)
    JMLR = Mul(JML,JMR)
    if (LMIJ /= JMLR) cycle
    W0TT1 = W0_TT(1)
    W1TT1 = W1_TT(1)
    NI = mod(LRJ-LRI,2)
    if (NI == 0) then
      W0TT1 = -W0TT1
      W1TT1 = -W1TT1
    end if
    !LIST = LIST3(LRI,LRJ,LRB)
    VLOP0 = W0*W0TT1
    VLOP1 = W1*W1TT1
    if (LINE == 26) then    !LRI,LRJ,LRA
      call COMP_LOOP_G(5,LRI,LRJ,LRS,LRA,VLOP0,VLOP1,WL0,LIST0,WL1,LIST1)
    end if
    if (LINE == 28) then    !LRI,LRJ,LRS,LRA
      call COMP_LOOP_G(7,LRI,LRJ,LRS,LRA,VLOP0,VLOP1,WL0,LIST0,WL1,LIST1)
    end if
    if (LINE == 29) then    !LRI,LRJ,LRS,LRA
      call COMP_LOOP_G(6,LRI,LRJ,LRS,LRA,VLOP0,VLOP1,WL0,LIST0,WL1,LIST1)
    end if
    !WL = (VLOP0-VLOP1)*vint_ci(LIST)-2*VLOP0*vint_ci(LIST+1)
    ! TT(11-1) (22)Ar(23)-Bl(32)-
    do LRK=NORB_FRZ+1,LRI-1
      LMK = LSM_INN(LRK)
      LMKI = Mul(LMK,LMI)
      LMKJ = Mul(LMK,LMJ)
      if ((LMKI == JML) .and. (LMKJ == JMR)) then
        IWDL = JUST(LRK,LRI)
        IWDR = JUST(LRK,LRJ)
        !call PRODAB(3,JPEL,IWDL,IWDR,JWL,JWR,WL,JPER)
        !call PRODAB_1(3,JPEL,IWDL,IWDR,JWL,JWR,WL2,JPER,LIST2,LIST3)
        call PRODAB_2(3,JPEL,IWDL,IWDR,JWL,JWR,WL0,JPER,LIST0)
        if (LIST1 /= 0) call PRODAB_2(3,JPEL,IWDL,IWDR,JWL,JWR,WL1,JPER,LIST1)

      end if
    end do
    ! TT(11-1) Ar(23)-Bl(32)-C"(22)-    ACT -C"-
    do LRK=LRJ+1,NORB_DZ
      LMK = LSM_INN(LRK)
      LMKI = Mul(LMK,LMI)
      LMKJ = Mul(LMK,LMJ)
      if ((LMKI == JML) .and. (LMKJ == JMR)) then
        IWDL = JUST(LRI,LRK)
        IWDR = JUST(LRJ,LRK)
        !call PRODAB(3,JPEL,IWDL,IWDR,JWL,JWR,WL,JPER)
        !call PRODAB_1(3,JPEL,IWDL,IWDR,JWL,JWR,WL2,JPER,LIST2,LIST3)
        call PRODAB_2(3,JPEL,IWDL,IWDR,JWL,JWR,WL0,JPER,LIST0)
        if (LIST1 /= 0) call PRODAB_2(3,JPEL,IWDL,IWDR,JWL,JWR,WL1,JPER,LIST1)

      end if
    end do
    ! TT(11-1) Ar(23)-C'(22)-Bl(32)-    ACT -C"-
    do LRK=LRI+1,LRJ-1
      LMK = LSM_INN(LRK)
      LMKI = Mul(LMK,LMI)
      LMKJ = Mul(LMK,LMJ)
      if ((LMKI == JML) .and. (LMKJ == JMR)) then
        IWDL = JUST(LRI,LRK)
        IWDR = JUST(LRK,LRJ)
        !call PRODAB(3,JPEL,IWDL,IWDR,JWL,JWR,-WL,JPER)
        !call PRODAB_1(3,JPEL,IWDL,IWDR,JWL,JWR,-WL2,JPER,LIST2,LIST3
        call PRODAB_2(3,JPEL,IWDL,IWDR,JWL,JWR,-WL0,JPER,LIST0)
        if (LIST1 /= 0) call PRODAB_2(3,JPEL,IWDL,IWDR,JWL,JWR,-WL1,JPER,LIST1)

      end if
    end do
  end do
end do

if (JPAD /= JPADL) return
if (JWL >= JWR) return

W0TT2 = W0_TT(2)
W1TT2 = W1_TT(2)
W0TT3 = W0_TT(3)

do LRI=NORB_FRZ+1,NORB_DZ-1
  LMI = LSM_INN(LRI)
  do LRJ=LRI+1,NORB_DZ
    LMJ = LSM_INN(LRJ)
    LMIJ = Mul(LMI,LMJ)
    if ((LMIJ /= JML) .or. (LMIJ /= JMR)) cycle
    if (JWL >= JWR) cycle
    ! TT(11-2) (22)Drl(22)-
    ! TT(11-2) Drl(22)-C"(22)-
    IWDL = JUST(LRI,LRJ)
    IWDR = IWDL
    VLOP0 = W0*W0TT2
    VLOP1 = W1*W1TT2
    if (LINE == 26) then    !LRJ,LRA
      call COMP_LOOP_G(9,LRJ,0,0,LRA,VLOP0,VLOP1,WL0,LIST0,WL1,LIST1)
    end if
    if (LINE == 28) then    !LRJ,LRS,LRA
      call COMP_LOOP_G(12,LRJ,0,LRS,LRA,VLOP0,VLOP1,WL0,LIST0,WL1,LIST1)
    end if
    if (LINE == 29) then    !LRJ,LRS,LRA
      call COMP_LOOP_G(11,LRJ,0,LRS,LRA,VLOP0,VLOP1,WL0,LIST0,WL1,LIST1)
    end if

    !call PRODAB_1(3,JPEL,IWDL,IWDR,JWL,JWR,WL2,JPER,LIST2,LIST3)
    call PRODAB_2(3,JPEL,IWDL,IWDR,JWL,JWR,WL0,JPER,LIST0)
    if (LIST1 /= 0) call PRODAB_2(3,JPEL,IWDL,IWDR,JWL,JWR,WL1,JPER,LIST1)

    if (LINE == 26) then    !LRI,LRA
      call COMP_LOOP_G(9,LRI,0,0,LRA,VLOP0,VLOP1,WL0,LIST0,WL1,LIST1)
    end if
    if (LINE == 28) then    !LRI,LRS,LRA
      call COMP_LOOP_G(12,LRI,0,LRS,LRA,VLOP0,VLOP1,WL0,LIST0,WL1,LIST1)
    end if
    if (LINE == 29) then    !LRI,LRS,LRA
      call COMP_LOOP_G(11,LRI,0,LRS,LRA,VLOP0,VLOP1,WL0,LIST0,WL1,LIST1)
    end if
    !call PRODAB_1(3,JPEL,IWDL,IWDR,JWL,JWR,WL2,JPER,LIST2,LIST3)
    call PRODAB_2(3,JPEL,IWDL,IWDR,JWL,JWR,WL0,JPER,LIST0)
    if (LIST1 /= 0) call PRODAB_2(3,JPEL,IWDL,IWDR,JWL,JWR,WL1,JPER,LIST1)

    !WL = WL+WLTMP
    !WL = (VLOP0-VLOP1)*(VOINT(LRB,LRI)+VOINT(LRB,LRJ))
    VLOP0 = W0*W0TT3
    VLOP1 = Zero
    do LRK=1,NORB_DZ
      if (LRK == LRI) cycle
      if (LRK == LRJ) cycle
      ! TT(11-3) Drl(33)-C"(22)-C"(22)-
      ! TT(11-3) (22)Drl(33)-C"(22)-
      ! TT(11-3) (22)(22)Drl(33)-
      if (LINE == 26) then    !LRK,LRA
        call COMP_LOOP_G(9,LRK,0,0,LRA,VLOP0,VLOP1,WL0,LIST0,WL1,LIST1)
      end if
      if (LINE == 28) then    !LRK,LRS,LRA
        call COMP_LOOP_G(12,LRK,0,LRS,LRA,VLOP0,VLOP1,WL0,LIST0,WL1,LIST1)
      end if
      if (LINE == 29) then    !LRK,LRS,LRA
        call COMP_LOOP_G(11,LRK,0,LRS,LRA,VLOP0,VLOP1,WL0,LIST0,WL1,LIST1)
      end if
      !WL = WL+WLTMP
      !call PRODAB_1(3,JPEL,IWDL,IWDR,JWL,JWR,WL2,JPER,LIST2,LIST3)
      call PRODAB_2(3,JPEL,IWDL,IWDR,JWL,JWR,WL0,JPER,LIST0)
      if (LIST1 /= 0) call PRODAB_2(3,JPEL,IWDL,IWDR,JWL,JWR,WL1,JPER,LIST1)

    end do
    !call PRODAB(3,JPEL,IWDL,IWDR,JWL,JWR,WL,JPER)
  end do
end do

return

end subroutine TT_HEAD_DBL_TAIL_ACT_G

subroutine TTTT_HEAD_DBL_TAIL_ACT_G(LRA)

use gugaci_global, only: jml, jmr, jpad, jpadl, jpel, jper, just, jwl, jwr, line, lrs, lsm_inn, norb_dz, norb_frz, w0, w0_t1t1, &
                         w1, w1_t1t1
use Symmetry_Info, only: Mul
use Constants, only: Zero
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: lra
integer(kind=iwp) :: iwdl, iwdr, list0, list1, lmi, lmij, lmj, lmm, lri, lrj, lrk, lrm, ni
real(kind=wp) :: vlop0, vlop1, w0tt1, w1tt1, wl0, wl1

! T1T1(12-1)  Ar(13)-Bl(31)-
! T1T1(12-1)  Ar(13)-C'(11)-Bl(31)-
! T1T1(12-1)  Ar(13)-Bl(31)-C"(11)-
do LRI=NORB_FRZ+1,NORB_DZ-1
  LMI = LSM_INN(LRI)
  do LRJ=LRI+1,NORB_DZ
    LMJ = LSM_INN(LRJ)
    W0TT1 = W0_T1T1(1)
    W1TT1 = W1_T1T1(1)
    NI = mod(LRJ-LRI,2)
    if (NI == 0) then
      W0TT1 = -W0TT1
      W1TT1 = -W1TT1
    end if
    !LIST = LIST4(LRI,LRJ,LRA,LRB)
    VLOP0 = W0*W0TT1
    VLOP1 = W1*W1TT1
    if (LINE == 26) then    !LRI,LRJ,LRA
      call COMP_LOOP_G(5,LRI,LRJ,LRS,LRA,VLOP0,VLOP1,WL0,LIST0,WL1,LIST1)
    end if
    if (LINE == 28) then    !LRI,LRJ,LRS,LRA
      call COMP_LOOP_G(7,LRI,LRJ,LRS,LRA,VLOP0,VLOP1,WL0,LIST0,WL1,LIST1)
    end if
    if (LINE == 29) then    !LRI,LRJ,LRS,LRA
      call COMP_LOOP_G(6,LRI,LRJ,LRS,LRA,VLOP0,VLOP1,WL0,LIST0,WL1,LIST1)
    end if
    !WL = (VLOP0-VLOP1)*VINT_CI(LIST)-2*VLOP0*VINT_CI(LIST+1)
    ! T1T1(12-1)  (11)Ar(13)-Bl(31)-
    do LRM=NORB_FRZ+1,LRI-1
      LMM = LSM_INN(LRM)
      if ((JML /= Mul(LMI,LMM)) .or. (JMR /= Mul(LMM,LMJ))) cycle
      IWDL = JUST(LRM,LRI)
      IWDR = JUST(LRM,LRJ)
      !call PRODAB(3,JPEL,IWDL,IWDR,JWL,JWR,WL,JPER)
      !call PRODAB_1(3,JPEL,IWDL,IWDR,JWL,JWR,WL2,JPER,LIST2,LIST3)
      call PRODAB_2(3,JPEL,IWDL,IWDR,JWL,JWR,WL0,JPER,LIST0)
      if (LIST1 /= 0) call PRODAB_2(3,JPEL,IWDL,IWDR,JWL,JWR,WL1,JPER,LIST1)

    end do
    ! T1T1(12-1)  Ar(13)-Bl(31)-C"(11)-
    do LRM=LRJ+1,NORB_DZ
      LMM = LSM_INN(LRM)
      if ((JML /= Mul(LMI,LMM)) .or. (JMR /= Mul(LMM,LMJ))) cycle
      IWDL = JUST(LRI,LRM)
      IWDR = JUST(LRJ,LRM)
      !call PRODAB(3,JPEL,IWDL,IWDR,JWL,JWR,WL,JPER)
      !call PRODAB_1(3,JPEL,IWDL,IWDR,JWL,JWR,WL2,JPER,LIST2,LIST3)
      call PRODAB_2(3,JPEL,IWDL,IWDR,JWL,JWR,WL0,JPER,LIST0)
      if (LIST1 /= 0) call PRODAB_2(3,JPEL,IWDL,IWDR,JWL,JWR,WL1,JPER,LIST1)

    end do
    ! T1T1(12-1)  Ar(13)-C'(11)-Bl(31)-
    do LRM=LRI+1,LRJ-1
      LMM = LSM_INN(LRM)
      if ((JML /= Mul(LMI,LMM)) .or. (JMR /= Mul(LMM,LMJ))) cycle
      IWDL = JUST(LRI,LRM)
      IWDR = JUST(LRM,LRJ)
      !call PRODAB(3,JPEL,IWDL,IWDR,JWL,JWR,-WL,JPER)
      !call PRODAB_1(3,JPEL,IWDL,IWDR,JWL,JWR,-WL2,JPER,LIST2,LIST3)
      call PRODAB_2(3,JPEL,IWDL,IWDR,JWL,JWR,-WL0,JPER,LIST0)
      if (LIST1 /= 0) call PRODAB_2(3,JPEL,IWDL,IWDR,JWL,JWR,-WL1,JPER,LIST1)

    end do
  end do
end do

if (JPAD /= JPADL) return
if (JWL >= JWR) return

do LRI=NORB_FRZ+1,NORB_DZ
  LMI = LSM_INN(LRI)
  do LRJ=LRI+1,NORB_DZ               !BBS_TMP
    LMJ = LSM_INN(LRJ)
    LMIJ = Mul(LMI,LMJ)
    if (JML /= LMIJ) cycle
    VLOP0 = W0*W0_T1T1(2)
    VLOP1 = W1*W1_T1T1(2)
    !-------------------------------------------------------------------
    ! lyb
    IWDL = JUST(LRI,LRJ)
    IWDR = IWDL

    ! T1T1(12-2)  (11)Drl(11)-
    if (LINE == 26) then    !LRI,LRA
      call COMP_LOOP_G(9,LRI,0,0,LRA,VLOP0,VLOP1,WL0,LIST0,WL1,LIST1)
    end if
    if (LINE == 28) then    !LRK,LRS,LRA
      call COMP_LOOP_G(12,LRI,0,LRS,LRA,VLOP0,VLOP1,WL0,LIST0,WL1,LIST1)
    end if
    if (LINE == 29) then    !LRK,LRS,LRA
      call COMP_LOOP_G(11,LRI,0,LRS,LRA,VLOP0,VLOP1,WL0,LIST0,WL1,LIST1)
    end if
    !call PRODAB_1(3,JPEL,IWDL,IWDR,JWL,JWR,WL2,JPER,LIST2,LIST3)
    call PRODAB_2(3,JPEL,IWDL,IWDR,JWL,JWR,WL0,JPER,LIST0)
    if (LIST1 /= 0) call PRODAB_2(3,JPEL,IWDL,IWDR,JWL,JWR,WL1,JPER,LIST1)

    ! T1T1(12-2)  Drl(11)-C"(11)-
    if (LINE == 26) then    !LRJ,LRA
      call COMP_LOOP_G(9,LRJ,0,0,LRA,VLOP0,VLOP1,WL0,LIST0,WL1,LIST1)
    end if
    if (LINE == 28) then    !LRJ,LRS,LRA
      call COMP_LOOP_G(12,LRJ,0,LRS,LRA,VLOP0,VLOP1,WL0,LIST0,WL1,LIST1)
    end if
    if (LINE == 29) then    !LRJ,LRS,LRA
      call COMP_LOOP_G(11,LRJ,0,LRS,LRA,VLOP0,VLOP1,WL0,LIST0,WL1,LIST1)
    end if
    !WL = WL+WLTMP
    !call PRODAB_1(3,JPEL,IWDL,IWDR,JWL,JWR,WL2,JPER,LIST2,LIST3)
    call PRODAB_2(3,JPEL,IWDL,IWDR,JWL,JWR,WL0,JPER,LIST0)
    if (LIST1 /= 0) call PRODAB_2(3,JPEL,IWDL,IWDR,JWL,JWR,WL1,JPER,LIST1)

    ! T1T1(12-3)  (11)(11)Drl(33)-
    ! T1T1(12-3)  (11)Drl(33)-C"(11)-
    ! T1T1(12-3)  Drl(33)-C"(11)-C"(11)-
    do LRK=1,NORB_DZ
      if (LRK == LRI) cycle
      if (LRK == LRJ) cycle
      VLOP0 = W0*W0_T1T1(3)
      VLOP1 = Zero
      if (LINE == 26) then    !LRK,LRA
        call COMP_LOOP_G(9,LRK,0,0,LRA,VLOP0,VLOP1,WL0,LIST0,WL1,LIST1)
      end if
      if (LINE == 28) then    !LRK,LRS,LRA
        call COMP_LOOP_G(12,LRK,0,LRS,LRA,VLOP0,VLOP1,WL0,LIST0,WL1,LIST1)
      end if
      if (LINE == 29) then    !LRK,LRS,LRA
        call COMP_LOOP_G(11,LRK,0,LRS,LRA,VLOP0,VLOP1,WL0,LIST0,WL1,LIST1)
      end if
      !WL = WL+WLTMP
      !call PRODAB_1(3,JPEL,IWDL,IWDR,JWL,JWR,WL2,JPER,LIST2,LIST3)
      call PRODAB_2(3,JPEL,IWDL,IWDR,JWL,JWR,WL0,JPER,LIST0)
      if (LIST1 /= 0) call PRODAB_2(3,JPEL,IWDL,IWDR,JWL,JWR,WL1,JPER,LIST1)

    end do
    !IWDL = JUST(LRI,LRJ)
    !IWDR = IWDL
    !call PRODAB(3,JPEL,IWDL,IWDR,JWL,JWR,WL,JPER)
  end do
end do

return

end subroutine TTTT_HEAD_DBL_TAIL_ACT_G

subroutine DD_HEAD_DBL_TAIL_ACT_G(LRA)

use gugaci_global, only: jml, jmr, jpel, jper, jud, jwl, jwr, line, lrs, lsm_inn, norb_dz, norb_frz, w0, w0_dd, w1, w1_dd
use Constants, only: Zero
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: lra
integer(kind=iwp) :: imil, imir, iwdl, iwdr, list0, list1, lril, lrir, lrk, ni
real(kind=wp) :: vlop0, vlop1, w0dd1, w1dd1, wl0, wl1

! DD(19-1) Ar(23)-Bl(32)-
! DD(19-2) Drl(22)-
! DD(19-3) Drl(33)-
! DD(19-3) Drl(33)-C"(22)-
do LRIL=NORB_FRZ+1,NORB_DZ
  IMIL = LSM_INN(LRIL)
  if (IMIL /= JML) cycle
  IWDL = JUD(LRIL)
  do LRIR=LRIL,NORB_DZ
    IMIR = LSM_INN(LRIR)
    if (IMIR /= JMR) cycle
    IWDR = JUD(LRIR)

    W0DD1 = W0_DD(1)
    W1DD1 = W1_DD(1)
    NI = mod(LRIR-LRIL,2)
    if (NI == 0) then
      W0DD1 = -W0DD1
      W1DD1 = -W1DD1
    end if

    if ((LRIL == LRIR) .and. (JWL < JWR)) then
      ! DD(19-2) Drl(22)-
      VLOP0 = W0*W0_DD(2)
      VLOP1 = W1*W1_DD(2)
      if (LINE == 26) then    !LRIL,LRA
        call COMP_LOOP_G(9,LRIL,0,0,LRA,VLOP0,VLOP1,WL0,LIST0,WL1,LIST1)
      end if
      if (LINE == 28) then    !LRIL,LRS,LRA
        call COMP_LOOP_G(12,LRIL,0,LRS,LRA,VLOP0,VLOP1,WL0,LIST0,WL1,LIST1)
      end if
      if (LINE == 29) then    !LRIL,LRS,LRA
        call COMP_LOOP_G(11,LRIL,0,LRS,LRA,VLOP0,VLOP1,WL0,LIST0,WL1,LIST1)
      end if
      !call PRODAB_1(3,JPEL,IWDL,IWDR,JWL,JWR,WL2,JPER,LIST2,LIST3)
      call PRODAB_2(3,JPEL,IWDL,IWDR,JWL,JWR,WL0,JPER,LIST0)
      if (LIST1 /= 0) call PRODAB_2(3,JPEL,IWDL,IWDR,JWL,JWR,WL1,JPER,LIST1)

      !WL = (VLOP0-VLOP1)*VOINT(LRA,LRIL)
      !call PRODAB(3,JPEL,IWDL,IWDR,JWL,JWR,WL,JPER)
      ! DD(19-3) Drl(33)-
      VLOP0 = W0*W0_DD(3)
      VLOP1 = Zero
      do LRK=1,NORB_DZ
        if (LRK == LRIL) cycle
        !LIST = LIST3(LRS,LRA,LRK)
        if (LINE == 26) then    !LRK,LRA
          call COMP_LOOP_G(9,LRK,0,0,LRA,VLOP0,VLOP1,WL0,LIST0,WL1,LIST1)     !WYB
        end if
        if (LINE == 28) then    !LRK,LRS,LRA
          call COMP_LOOP_G(12,LRK,0,LRS,LRA,VLOP0,VLOP1,WL0,LIST0,WL1,LIST1)   !WYB
        end if
        if (LINE == 29) then    !LRK,LRS,LRA
          call COMP_LOOP_G(11,LRK,0,LRS,LRA,VLOP0,VLOP1,WL0,LIST0,WL1,LIST1)   !WYB
        end if
        !call PRODAB_1(3,JPEL,IWDL,IWDR,JWL,JWR,WL2,JPER,LIST2,LIST3)
        call PRODAB_2(3,JPEL,IWDL,IWDR,JWL,JWR,WL0,JPER,LIST0)
        if (LIST1 /= 0) call PRODAB_2(3,JPEL,IWDL,IWDR,JWL,JWR,WL1,JPER,LIST1)
        !WL = WL+VLOP0*VOINT(LRK,LRA)
        !WL = WLTMP+WL
      end do
      !call PRODAB(3,JPEL,IWDL,IWDR,JWL,JWR,WL,JPER)
    end if
    if (LRIL /= LRIR) then
      ! DD(19-1) Ar(23)-Bl(32)-
      VLOP0 = W0*W0DD1
      VLOP1 = W1*W1DD1
      if (LINE == 26) then   !LRIL,LRIR,LRA
        call COMP_LOOP_G(5,LRIL,LRIR,LRS,LRA,VLOP0,VLOP1,WL0,LIST0,WL1,LIST1)
      end if
      if (LINE == 28) then   !LRIL,LRIR,LRS,LRA
        call COMP_LOOP_G(7,LRIL,LRIR,LRS,LRA,VLOP0,VLOP1,WL0,LIST0,WL1,LIST1)
      end if
      if (LINE == 29) then   !LRIL,LRIR,LRS,LRA
        call COMP_LOOP_G(6,LRIL,LRIR,LRS,LRA,VLOP0,VLOP1,WL0,LIST0,WL1,LIST1)
      end if
      !LIST = LIST4(LRIL,LRIR,LRS,LRA)
      !WL = VLOP0*(vint_ci(LIST)-2*vint_ci(LIST+1))-VLOP1*vint_ci(LIST) !1.1       !!!!!
      !call PRODAB(3,JPEL,IWDL,IWDR,JWL,JWR,WL,JPER)
      !call PRODAB_1(3,JPEL,IWDL,IWDR,JWL,JWR,WL2,JPER,LIST2,LIST3)
      call PRODAB_2(3,JPEL,IWDL,IWDR,JWL,JWR,WL0,JPER,LIST0)
      if (LIST1 /= 0) call PRODAB_2(3,JPEL,IWDL,IWDR,JWL,JWR,WL1,JPER,LIST1)

    end if
  end do
end do

return

end subroutine DD_HEAD_DBL_TAIL_ACT_G

subroutine DDDD_HEAD_DBL_TAIL_ACT_G(LRA)

use gugaci_global, only: jml, jmr, jpel, jper, jud, jwl, jwr, line, lrs, lsm_inn, norb_dz, norb_frz, w0, w0_d1d1, w1, w1_d1d1
use Constants, only: Zero
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: lra
integer(kind=iwp) :: imil, imir, iwdl, iwdr, list0, list1, lril, lrir, lrk, ni
real(kind=wp) :: vlop0, vlop1, w0dd1, w1dd1, wl0, wl1

! D1D1(20-1) Ar(13)-BL(31)-
! D1D1(20-1) Drl(11)-
! D1D1(20-1) Drl(33)-
! D1D1(20-1) Drl(33)-C"(11)-
do LRIL=NORB_FRZ+1,NORB_DZ
  IMIL = LSM_INN(LRIL)
  if (IMIL /= JML) cycle
  IWDL = JUD(LRIL)
  do LRIR=LRIL,NORB_DZ
    IMIR = LSM_INN(LRIR)
    if (IMIR /= JMR) cycle
    W0DD1 = W0_D1D1(1)
    W1DD1 = W1_D1D1(1)
    NI = mod(LRIR-LRIL,2)
    if (NI == 0) then
      W0DD1 = -W0DD1
      W1DD1 = -W1DD1
    end if
    IWDR = JUD(LRIR)
    if ((LRIL == LRIR) .and. (JWL < JWR)) then
      ! D1D1(20-1) Drl(11)-
      VLOP0 = W0*W0_D1D1(2)
      VLOP1 = W1*W1_D1D1(2)
      if (LINE == 26) then    !LRIL,LRA
        call COMP_LOOP_G(9,LRIL,0,0,LRA,VLOP0,VLOP1,WL0,LIST0,WL1,LIST1)
      end if
      if (LINE == 28) then    !LRIL,LRS,LRA
        call COMP_LOOP_G(12,LRIL,0,LRS,LRA,VLOP0,VLOP1,WL0,LIST0,WL1,LIST1)
      end if
      if (LINE == 29) then    !LRIL,LRS,LRA
        call COMP_LOOP_G(11,LRIL,0,LRS,LRA,VLOP0,VLOP1,WL0,LIST0,WL1,LIST1)
      end if
      !call PRODAB_1(3,JPEL,IWDL,IWDR,JWL,JWR,WL2,JPER,LIST2,LIST3)
      call PRODAB_2(3,JPEL,IWDL,IWDR,JWL,JWR,WL0,JPER,LIST0)
      if (LIST1 /= 0) call PRODAB_2(3,JPEL,IWDL,IWDR,JWL,JWR,WL1,JPER,LIST1)

      ! D1D1(20-1) Drl(33)-
      ! D1D1(20-1) Drl(33)-C"(11)-
      VLOP0 = W0*W0_D1D1(3)
      VLOP1 = Zero
      do LRK=1,NORB_DZ
        if (LRK == LRIL) cycle
        if (LINE == 26) then    !LRK,LRA
          call COMP_LOOP_G(9,LRK,0,0,LRA,VLOP0,VLOP1,WL0,LIST0,WL1,LIST1)
        end if
        if (LINE == 28) then    !LRK,LRS,LRA
          call COMP_LOOP_G(12,LRK,0,LRS,LRA,VLOP0,VLOP1,WL0,LIST0,WL1,LIST1)
        end if
        if (LINE == 29) then    !LRK,LRS,LRA
          call COMP_LOOP_G(11,LRK,0,LRS,LRA,VLOP0,VLOP1,WL0,LIST0,WL1,LIST1)
        end if
        !WL = WL+WLTMP
        !call PRODAB_1(3,JPEL,IWDL,IWDR,JWL,JWR,WL2,JPER,LIST2,LIST3)
        call PRODAB_2(3,JPEL,IWDL,IWDR,JWL,JWR,WL0,JPER,LIST0)
        if (LIST1 /= 0) call PRODAB_2(3,JPEL,IWDL,IWDR,JWL,JWR,WL1,JPER,LIST1)

      end do
      !call PRODAB(3,JPEL,IWDL,IWDR,JWL,JWR,WL,JPER)
    end if
    if (LRIL /= LRIR) then
      ! D1D1(20-1) Ar(13)-BL(31)-
      VLOP0 = W0*W0DD1
      VLOP1 = W1*W1DD1
      if (LINE == 26) then   !LRIL,LRIR,LRA
        call COMP_LOOP_G(5,LRIL,LRIR,LRS,LRA,VLOP0,VLOP1,WL0,LIST0,WL1,LIST1)
      end if
      if (LINE == 28) then   !LRIL,LRIR,LRS,LRA
        call COMP_LOOP_G(7,LRIL,LRIR,LRS,LRA,VLOP0,VLOP1,WL0,LIST0,WL1,LIST1)
      end if
      if (LINE == 29) then   !LRIL,LRIR,LRS,LRA
        call COMP_LOOP_G(6,LRIL,LRIR,LRS,LRA,VLOP0,VLOP1,WL0,LIST0,WL1,LIST1)
      end if
      !LIST = LIST3(LRIL,LRIR,LRA)
      !WL = (VLOP0-VLOP1)*vint_ci(LIST)-2*VLOP0*vint_ci(LIST+1)   !2
      !call PRODAB(3,JPEL,IWDL,IWDR,JWL,JWR,WL,JPER)
      !call PRODAB_1(3,JPEL,IWDL,IWDR,JWL,JWR,WL2,JPER,LIST2,LIST3)
      call PRODAB_2(3,JPEL,IWDL,IWDR,JWL,JWR,WL0,JPER,LIST0)
      if (LIST1 /= 0) call PRODAB_2(3,JPEL,IWDL,IWDR,JWL,JWR,WL1,JPER,LIST1)

    end if
  end do
end do

return

end subroutine DDDD_HEAD_DBL_TAIL_ACT_G

subroutine DD1_HEAD_DBL_TAIL_ACT_G(LRA)

use gugaci_global, only: jml, jmr, jpel, jper, jud, jwl, jwr, line, lrs, lsm_inn, norb_dz, norb_frz, w1, w1_dd1
use Constants, only: Zero
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: lra
integer(kind=iwp) :: iwdl, iwdr, list0, list1, lmi, lmj, lri, lrj, ni
real(kind=wp) :: vlop0, vlop1, wl0, wl1

do LRI=NORB_FRZ+1,NORB_DZ
  LMI = LSM_INN(LRI)
  IWDL = JUD(LRI)
  do LRJ=LRI+1,NORB_DZ
    ! DD1(21) Ar(23)-Bl(31)-
    LMJ = LSM_INN(LRJ)
    if ((JML /= LMI) .or. (JMR /= LMJ)) cycle
    VLOP1 = W1*W1_DD1
    VLOP0 = Zero
    NI = mod(LRJ-LRI,2)
    if (NI == 0) then
      VLOP1 = -VLOP1
    end if

    if (LINE == 26) then   !LRI,LRJ,LRA
      call COMP_LOOP_G(5,LRI,LRJ,LRS,LRA,VLOP0,VLOP1,WL0,LIST0,WL1,LIST1)
    end if
    if (LINE == 28) then   !LRIL,LRJ,LRS,LRA
      call COMP_LOOP_G(7,LRI,LRJ,LRS,LRA,VLOP0,VLOP1,WL0,LIST0,WL1,LIST1)
    end if
    if (LINE == 29) then   !LRI,LRJ,LRS,LRA
      call COMP_LOOP_G(6,LRI,LRJ,LRS,LRA,VLOP0,VLOP1,WL0,LIST0,WL1,LIST1)
    end if
    !LIST = LIST3(LRI,LRJ,LRA)
    !WL = -VLOP1*VINT_CI(LIST)
    IWDR = JUD(LRJ)
    !call PRODAB(3,JPEL,IWDL,IWDR,JWL,JWR,WL,JPER)
    !call PRODAB_1(3,JPEL,IWDL,IWDR,JWL,JWR,WL2,JPER,LIST2,LIST3)
    call PRODAB_2(3,JPEL,IWDL,IWDR,JWL,JWR,WL0,JPER,LIST0)
    if (LIST1 /= 0) call PRODAB_2(3,JPEL,IWDL,IWDR,JWL,JWR,WL1,JPER,LIST1)

  end do
end do

return

end subroutine DD1_HEAD_DBL_TAIL_ACT_G

subroutine D1D_HEAD_DBL_TAIL_ACT_G(LRA)

use gugaci_global, only: jml, jmr, jpel, jper, jud, jwl, jwr, line, lrs, lsm_inn, norb_dz, norb_frz, w1, w1_d1d
use Constants, only: Zero
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: lra
integer(kind=iwp) :: iwdl, iwdr, list0, list1, lmi, lmj, lri, lrj
real(kind=wp) :: vlop0, vlop1, wl0, wl1

! D1D(22-1)   Ar(13)-Bl(32)-
! D1D(22-2)   Drl(12)-
do LRI=NORB_FRZ+1,NORB_DZ-1
  LMI = LSM_INN(LRI)
  IWDL = JUD(LRI)
  do LRJ=LRI+1,NORB_DZ
    LMJ = LSM_INN(LRJ)
    if ((JML /= LMI) .or. (JMR /= LMJ)) cycle
    VLOP1 = W1*W1_D1D(1)
    VLOP0 = Zero
    if (mod(LRJ-LRI,2) == 0) then
      VLOP1 = -VLOP1
    end if
    if (LINE == 26) then   !LRI,LRJ,LRA
      call COMP_LOOP_G(5,LRI,LRJ,LRS,LRA,VLOP0,VLOP1,WL0,LIST0,WL1,LIST1)
    end if
    if (LINE == 28) then   !LRIL,LRJ,LRS,LRA
      call COMP_LOOP_G(7,LRI,LRJ,LRS,LRA,VLOP0,VLOP1,WL0,LIST0,WL1,LIST1)
    end if
    if (LINE == 29) then   !LRI,LRJ,LRS,LRA
      call COMP_LOOP_G(6,LRI,LRJ,LRS,LRA,VLOP0,VLOP1,WL0,LIST0,WL1,LIST1)
    end if
    !LIST = LIST3(LRI,LRJ,LRA)
    !WL = -VLOP1*VINT_CI(LIST)
    IWDR = JUD(LRJ)
    !call PRODAB(3,JPEL,IWDL,IWDR,JWL,JWR,WL,JPER)
    !call PRODAB_1(3,JPEL,IWDL,IWDR,JWL,JWR,WL2,JPER,LIST2,LIST3)
    call PRODAB_2(3,JPEL,IWDL,IWDR,JWL,JWR,WL0,JPER,LIST0)
    if (LIST1 /= 0) call PRODAB_2(3,JPEL,IWDL,IWDR,JWL,JWR,WL1,JPER,LIST1)

  end do
end do
VLOP1 = W1*W1_D1D(2)
VLOP0 = Zero
do LRI=NORB_FRZ+1,NORB_DZ
  ! D1D(22-2)   Drl(12)-
  LMI = LSM_INN(LRI)
  if ((JML /= LMI) .or. (JMR /= LMI)) cycle
  if (LINE == 26) then    !LRI,LRA
    call COMP_LOOP_G(9,LRI,0,0,LRA,VLOP0,VLOP1,WL0,LIST0,WL1,LIST1)
  end if
  if (LINE == 28) then    !LRI,LRS,LRA
    call COMP_LOOP_G(12,LRI,0,LRS,LRA,VLOP0,VLOP1,WL0,LIST0,WL1,LIST1)
  end if
  if (LINE == 29) then    !LRI,LRS,LRA
    call COMP_LOOP_G(11,LRI,0,LRS,LRA,VLOP0,VLOP1,WL0,LIST0,WL1,LIST1)
  end if
  !VLOP1 = W1*W1_D1D(2)
  !WL = -VLOP1*VOINT(LRI,LRA)
  IWDL = JUD(LRI)
  IWDR = IWDL
  !call PRODAB(3,JPEL,IWDL,IWDR,JWL,JWR,WL,JPER)
  !call PRODAB_1(3,JPEL,IWDL,IWDR,JWL,JWR,WL2,JPER,LIST2,LIST3)
  call PRODAB_2(3,JPEL,IWDL,IWDR,JWL,JWR,WL0,JPER,LIST0)
  if (LIST1 /= 0) call PRODAB_2(3,JPEL,IWDL,IWDR,JWL,JWR,WL1,JPER,LIST1)

end do

return

end subroutine D1D_HEAD_DBL_TAIL_ACT_G

subroutine SV_HEAD_DBL_TAIL_ACT_G(LRA)

use gugaci_global, only: jb_sys, jml, jpel, jper, just, jwl, jwr, line, lrs, lsm_inn, norb_dz, norb_frz, w0, w0_sv, w1, w1_sv
use Symmetry_Info, only: Mul
use Constants, only: Zero
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: lra
integer(kind=iwp) :: iwdl, iwdr, list0, list1, lmi, lmij, lmj, lri, lrj, ni
real(kind=wp) :: vlop0, vlop1, w0sv1, w0sv2, w1sv1, w1sv2, wl0, wl1

! SV(10-1) Ar(13)-Br(23)-
! SV(10-2) Ar(23)-Br(13)-
! SV(10-3) Drr(03)-
IWDR = 0
do LRI=NORB_FRZ+1,NORB_DZ
  LMI = LSM_INN(LRI)
  do LRJ=LRI,NORB_DZ
    LMJ = LSM_INN(LRJ)
    LMIJ = Mul(LMI,LMJ)
    if (LMIJ /= JML) cycle
    W0SV1 = W0_SV(1)
    W1SV1 = W1_SV(1)
    W0SV2 = W0_SV(2)
    W1SV2 = W1_SV(2)
    NI = mod(LRJ-LRI,2)
    if (NI == 0) then
      W0SV1 = -W0SV1
      W1SV1 = -W1SV1
      W0SV2 = -W0SV2
      W1SV2 = -W1SV2
    end if
    IWDL = JUST(LRI,LRJ)
    if (LRI == LRJ) then
      VLOP0 = W0*W0_SV(3)            !D10-3
      VLOP1 = Zero
      !WL = VLOP0*VOINT(LRA,LRI)/2
      if (LINE == 25) then    !LRI,LRA
        ! wyb call COMP_LOOP_G(8,LRI,LRG,LRS,LRA,VLOP0,VLOP1,WL)
        call COMP_LOOP_G(8,LRI,0,0,LRA,VLOP0,VLOP1,WL0,LIST0,WL1,LIST1)
      end if
      if (LINE == 27) then    !LRI,LRS,LRA
        ! wyb call COMP_LOOP_G(10,LRI,LRS,LRA,LRA,VLOP0,VLOP1,WL)
        call COMP_LOOP_G(10,LRI,0,LRS,LRA,VLOP0,VLOP1,WL0,LIST0,WL1,LIST1)
      end if
      !call PRODAB(3,JPEL,IWDL,IWDR,JWL,JWR,WL,JPER)
      !call PRODAB_1(3,JPEL,IWDL,IWDR,JWL,JWR,WL2,JPER,LIST2,LIST3)
      call PRODAB_2(3,JPEL,IWDL,IWDR,JWL,JWR,WL0,JPER,LIST0)
      if (LIST1 /= 0) call PRODAB_2(3,JPEL,IWDL,IWDR,JWL,JWR,WL1,JPER,LIST1)

    else
      VLOP0 = W0*W0SV2             !D10-2
      VLOP1 = W1*W1SV2
      !LIST = LIST3(LRI,LRJ,LRA)
      !WL = (VLOP0+VLOP1)*vint_ci(LIST)        !2.1          !!!!!
      if (LINE == 25) then    !LRI,LRJ,LRA
        call COMP_LOOP_G(3,LRI,LRJ,LRS,LRA,VLOP0,VLOP1,WL0,LIST0,WL1,LIST1)
      end if
      if (LINE == 27) then    !LRI,LRJ,LRS,LRA
        call COMP_LOOP_G(4,LRI,LRJ,LRS,LRA,VLOP0,VLOP1,WL0,LIST0,WL1,LIST1)
      end if
      !call PRODAB(3,JPEL,IWDL,IWDR,JWL,JWR,WL,JPER)
      !call PRODAB_1(3,JPEL,IWDL,IWDR,JWL,JWR,WL2,JPER,LIST2,LIST3)
      call PRODAB_2(3,JPEL,IWDL,IWDR,JWL,JWR,WL0,JPER,LIST0)
      if (LIST1 /= 0) call PRODAB_2(3,JPEL,IWDL,IWDR,JWL,JWR,WL1,JPER,LIST1)

      if (JB_SYS > 0) then
        IWDL = JUST(LRJ,LRI)
        VLOP0 = W0*W0SV1
        VLOP1 = W1*W1SV1
        !LIST = LIST3(LRI,LRJ,LRA)
        !WL = (VLOP0+VLOP1)*vint_ci(LIST)        !2.1          !!!!!
        if (LINE == 25) then    !LRI,LRJ,LRA
          call COMP_LOOP_G(3,LRI,LRJ,LRS,LRA,VLOP0,VLOP1,WL0,LIST0,WL1,LIST1)
        end if
        if (LINE == 27) then    !LRI,LRJ,LRS,LRA
          call COMP_LOOP_G(4,LRI,LRJ,LRS,LRA,VLOP0,VLOP1,WL0,LIST0,WL1,LIST1)
        end if
        !call PRODAB(3,JPEL,IWDL,IWDR,JWL,JWR,WL,JPER)
        !call PRODAB_1(3,JPEL,IWDL,IWDR,JWL,JWR,WL2,JPER,LIST2,LIST3)
        call PRODAB_2(3,JPEL,IWDL,IWDR,JWL,JWR,WL0,JPER,LIST0)
        if (LIST1 /= 0) call PRODAB_2(3,JPEL,IWDL,IWDR,JWL,JWR,WL1,JPER,LIST1)

      end if
    end if
  end do
end do

return

end subroutine SV_HEAD_DBL_TAIL_ACT_G

subroutine SD_HEAD_DBL_TAIL_ACT_G(LRA,LPCOE)
!**********************************************
!     LRA, ....partial loop.......
!**********************************************

use gugaci_global, only: jb_sys, jml, jmr, jpel, jper, jud, just, jwl, jwr, lsm_inn, norb_dz, norb_frz, norb_inn, w0, w0_sd, w1, &
                         w1_sd
use Symmetry_Info, only: Mul
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: lra, lpcoe(norb_dz+1:norb_inn)
integer(kind=iwp) :: iwdl, iwdl1, iwdr, jmlr, kcoe, lmd, lmi, lmij, lmj, lr, lrd, lri, lrj, lrk, ni, nocc, nxo
real(kind=wp) :: tcoe, vlop0, vlop1, w0sd1, w0sd11, w0sd12, w0sd14, w0sd15, w0sd16, w0sd2, w0sd3, w0sd4, w0sd5, w0sd6, w0sd7, &
                 w0sd8, w0sd9, w1sd10, w1sd11, w1sd2, w1sd5, w1sd6, w1sd7, w1sd8, wl

JMLR = Mul(JML,JMR)
! SD(6-1) A&r(02)-
! SD(6-2) C(22)A&r(13)-
! SD(6-4) A&r(23)C'(12)-
! SD(6-5) A&r(23)B&r(13)B^r(32)
! SD(6-8) A&r(23)B&l(32)B^l(13)
! SD(6-9) D&r&r(03)B^r(32)
! SD(6-11) D&r&l(22)B^l(13)
! SD(6-12) D&r&l(33)B^l(02)
! SD(6-13) (22)D&r&l(33)B^l(13)
! SD(6-14) D&r&l(33)C"(22)B^l(13)
! SD(6-16) D&r&l(33)B^l(23)C'(12)

! SD(6-3) A&r(13)C'(22)-
! SD(6-6) A&r(13)B&r(23)B^r(32)
! SD(6-7) A&r(13)B&l(32)B^l(23)
! SD(6-10) D&r&l(12)B^l(23)
! SD(6-15) D&r&l(33)B^l(13)C'(22)

if (JML == 1) then

  !SD(6-1) A&r(02)-
  do LRI=NORB_FRZ+1,NORB_DZ
    LMI = LSM_INN(LRI)
    IWDL = JUST(LRI,LRI)
    W0SD1 = W0_SD(1)
    W0SD12 = W0_SD(12)
    NI = mod(NORB_DZ-LRI,2)
    if (NI == 1) then
      W0SD1 = -W0SD1
      W0SD12 = -W0SD12
    end if
    if (LMI == JMLR) then
      IWDR = JUD(LRI)
      VLOP0 = W0*W0SD1
      !WL = VLOP0*VOINT(LRI,LRA)
      WL = VLOP0
      call PRODAB_1(3,JPEL,IWDL,IWDR,JWL,JWR,WL,JPER,LRI,LRA)

      do LR=LRI+1,NORB_DZ
        !LIST = LIST3(LRI,LRA,LR)
        !WL = WL+VLOP0*(2*VINT_CI(LIST+1)-VINT_CI(LIST))  ! 310:NEOC=2
        WL = VLOP0*2
        call TRANS_IJKL_INTPOS(LRA,LRI,LR,LR,NXO)
        call PRODAB_2(3,JPEL,IWDL,IWDR,JWL,JWR,WL,JPER,NXO)
        WL = -VLOP0
        call TRANS_IJKL_INTPOS(LRA,LR,LRI,LR,NXO)
        call PRODAB_2(3,JPEL,IWDL,IWDR,JWL,JWR,WL,JPER,NXO)
      end do
      do LRK=norb_dz+1,LRA
        !LIST = LIST3(LRI,LRA,LRK)
        KCOE = LPCOE(LRK)
        call NEOC(KCOE,NOCC,TCOE)
        !WL = WL+VLOP0*NOCC*(VINT_CI(LIST+1)+TCOE*VINT_CI(LIST))
        WL = VLOP0*NOCC
        call TRANS_IJKL_INTPOS(LRA,LRI,LRK,LRK,NXO)
        call PRODAB_2(3,JPEL,IWDL,IWDR,JWL,JWR,WL,JPER,NXO)
        WL = VLOP0*NOCC*TCOE
        call TRANS_IJKL_INTPOS(LRA,LRK,LRI,LRK,NXO)
        call PRODAB_2(3,JPEL,IWDL,IWDR,JWL,JWR,WL,JPER,NXO)

      end do
      !WL = WL*VLOP0
      ! SD(6-12) D&rl(33)B^l(02)
      VLOP0 = W0*W0SD12
      do LRK=1,LRI-1
        !LIST = LIST3(LRI,LRA,LRK)
        !WL = WL+VLOP0*(vint_ci(LIST)-2*VINT_CI(LIST+1))
        WL = VLOP0
        call TRANS_IJKL_INTPOS(LRA,LRK,LRI,LRK,NXO)
        call PRODAB_2(3,JPEL,IWDL,IWDR,JWL,JWR,WL,JPER,NXO)
        WL = -VLOP0*2
        call TRANS_IJKL_INTPOS(LRA,LRI,LRK,LRK,NXO)
        call PRODAB_2(3,JPEL,IWDL,IWDR,JWL,JWR,WL,JPER,NXO)

      end do
      !call PRODAB(3,JPEL,IWDL,IWDR,JWL,JWR,WL,JPER)
    end if

    ! SD(6-9) D&r&r(03)B^r(32)
    do LRD=LRI+1,NORB_DZ
      LMD = LSM_INN(LRD)
      if (LMD /= JMR) cycle
      W0SD9 = W0_SD(9)
      NI = mod(NORB_DZ-LRD,2)
      if (NI == 1) W0SD9 = -W0SD9
      IWDR = JUD(LRD)
      VLOP0 = W0*W0SD9
      !LIST = LIST3(LRD,LRA,LRI)
      !WL = VINT_CI(LIST)*VLOP0
      WL = VLOP0
      call TRANS_IJKL_INTPOS(LRA,LRI,LRD,LRI,NXO)
      call PRODAB_2(3,JPEL,IWDL,IWDR,JWL,JWR,WL,JPER,NXO)
    end do
  end do

end if

do LRI=NORB_FRZ+1,NORB_DZ
  LMI = LSM_INN(LRI)
  do LRJ=LRI+1,NORB_DZ
    LMJ = LSM_INN(LRJ)
    LMIJ = Mul(LMI,LMJ)
    if (LMIJ /= JML) cycle
    W0SD2 = W0_SD(2)
    W1SD2 = W1_SD(2)
    W0SD11 = W0_SD(11)
    W1SD11 = W1_SD(11)
    W0SD14 = W0_SD(14)
    NI = mod(NORB_DZ-LRJ,2)
    if (NI == 1) then
      W0SD2 = -W0SD2
      W1SD2 = -W1SD2
      W0SD11 = -W0SD11
      W1SD11 = -W1SD11
      W0SD14 = -W0SD14
    end if
    IWDL = JUST(LRI,LRJ)
    IWDL1 = JUST(LRJ,LRI)
    !*******************************************************************
    ! SD(6-2) C(22)-A&r(13)-
    if (LMI == JMR) then
      IWDR = JUD(LRI)
      VLOP0 = W0*W0SD2
      !LIST = LIST3(LRJ,LRA,LRJ)
      !WL = VLOP0*(VOINT(LRJ,LRA)+VINT_CI(LIST))          !310,act_c
      WL = VLOP0
      call TRANS_IJKL_INTPOS(LRA,LRJ,LRJ,LRJ,NXO)
      call PRODAB_1(3,JPEL,IWDL,IWDR,JWL,JWR,WL,JPER,LRJ,LRA)
      call PRODAB_2(3,JPEL,IWDL,IWDR,JWL,JWR,WL,JPER,NXO)

      do LR=LRJ+1,NORB_DZ
        !LIST = LIST3(LRJ,LRA,LR)
        !WL = WL+VLOP0*(2*VINT_CI(LIST+1)-VINT_CI(LIST)) !  310:NEOC
        WL = VLOP0*2
        call TRANS_IJKL_INTPOS(LRA,LRJ,LR,LR,NXO)
        call PRODAB_2(3,JPEL,IWDL,IWDR,JWL,JWR,WL,JPER,NXO)
        WL = -VLOP0
        call TRANS_IJKL_INTPOS(LRA,LR,LRJ,LR,NXO)
        call PRODAB_2(3,JPEL,IWDL,IWDR,JWL,JWR,WL,JPER,NXO)

      end do
      do LRK=norb_dz+1,LRA
        !LIST = LIST3(LRJ,LRA,LRK)
        KCOE = LPCOE(LRK)
        call NEOC(KCOE,NOCC,TCOE)
        !WL = WL+VLOP0*NOCC*(VINT_CI(LIST+1)+TCOE*VINT_CI(LIST))
        WL = VLOP0*NOCC
        call TRANS_IJKL_INTPOS(LRA,LRJ,LRK,LRK,NXO)
        call PRODAB_2(3,JPEL,IWDL,IWDR,JWL,JWR,WL,JPER,NXO)
        WL = VLOP0*NOCC*TCOE
        call TRANS_IJKL_INTPOS(LRA,LRK,LRJ,LRK,NXO)
        call PRODAB_2(3,JPEL,IWDL,IWDR,JWL,JWR,WL,JPER,NXO)

      end do
      !WL = WL*VLOP0
      ! SD(6-11) D&r&l(22)B^l(13)
      VLOP0 = W0*W0SD11
      VLOP1 = W1*W1SD11
      !LIST = LIST3(LRJ,LRA,LRI)
      !WL = WL+(VLOP0-VLOP1)*vint_ci(LIST)-2*VLOP0*VINT_CI(LIST+1)
      WL = VLOP0-VLOP1
      call TRANS_IJKL_INTPOS(LRA,LRI,LRJ,LRI,NXO)
      call PRODAB_2(3,JPEL,IWDL,IWDR,JWL,JWR,WL,JPER,NXO)
      WL = -VLOP0*2
      call TRANS_IJKL_INTPOS(LRA,LRJ,LRI,LRI,NXO)
      call PRODAB_2(3,JPEL,IWDL,IWDR,JWL,JWR,WL,JPER,NXO)

      ! SD(6-13) (22)D&r&l(33)B^l(13)
      ! SD(6-14) D&r&l(33)C"(22)B^l(13)
      VLOP0 = W0*W0SD14
      do LRK=1,LRJ-1
        if (LRK == LRI) cycle
        !LIST = LIST3(LRJ,LRA,LRK)
        !WL = WL+VLOP0*(vint_ci(LIST)-2*VINT_CI(LIST+1))
        WL = VLOP0
        call TRANS_IJKL_INTPOS(LRA,LRK,LRJ,LRK,NXO)
        call PRODAB_2(3,JPEL,IWDL,IWDR,JWL,JWR,WL,JPER,NXO)
        WL = -VLOP0*2
        call TRANS_IJKL_INTPOS(LRA,LRJ,LRK,LRK,NXO)
        call PRODAB_2(3,JPEL,IWDL,IWDR,JWL,JWR,WL,JPER,NXO)

      end do
      !call PRODAB(3,JPEL,IWDL,IWDR,JWL,JWR,WL,JPER)
    end if
    if (JB_SYS > 0) then
      if (LMJ == JMR) then
        IWDR = JUD(LRJ)
        W0SD3 = W0_SD(3)
        W0SD15 = W0_SD(15)
        W1SD10 = W1_SD(10)
        NI = mod(NORB_DZ-LRI,2)
        if (NI == 0) then
          W0SD3 = -W0SD3
          W0SD15 = -W0SD15
          W1SD10 = -W1SD10
        end if
        ! SD(6-3) A&r(13)-C'(22)-
        VLOP0 = W0*W0SD3
        !LIST = LIST3(LRI,LRA,LRI)
        !WL = VOINT(LRI,LRA)+VINT_CI(LIST)
        WL = VLOP0
        call PRODAB_1(3,JPEL,IWDL1,IWDR,JWL,JWR,WL,JPER,LRI,LRA)
        call TRANS_IJKL_INTPOS(LRA,LRI,LRI,LRI,NXO)
        call PRODAB_2(3,JPEL,IWDL1,IWDR,JWL,JWR,WL,JPER,NXO)
        !LIST = LIST3(LRI,LRA,LRJ)
        !WL = WL+VINT_CI(LIST+1)-VINT_CI(LIST)   !310 C'(22)NEOC=1,COE
        WL = VLOP0
        call TRANS_IJKL_INTPOS(LRA,LRI,LRJ,LRJ,NXO)
        call PRODAB_2(3,JPEL,IWDL1,IWDR,JWL,JWR,WL,JPER,NXO)
        WL = -VLOP0
        call TRANS_IJKL_INTPOS(LRA,LRJ,LRI,LRJ,NXO)
        call PRODAB_2(3,JPEL,IWDL1,IWDR,JWL,JWR,WL,JPER,NXO)

        do LR=LRI+1,NORB_DZ
          if (LR == LRJ) cycle
          !LIST = LIST3(LRI,LRA,LR)
          !WL = WL+2*VINT_CI(LIST+1)-VINT_CI(LIST)       !310:NEOC=2,C
          WL = VLOP0*2
          call TRANS_IJKL_INTPOS(LRA,LRI,LR,LR,NXO)
          call PRODAB_2(3,JPEL,IWDL1,IWDR,JWL,JWR,WL,JPER,NXO)
          WL = -VLOP0
          call TRANS_IJKL_INTPOS(LRA,LR,LRI,LR,NXO)
          call PRODAB_2(3,JPEL,IWDL1,IWDR,JWL,JWR,WL,JPER,NXO)

        end do
        do LRK=norb_dz+1,LRA
          !LIST = LIST3(LRI,LRA,LRK)
          KCOE = LPCOE(LRK)
          call NEOC(KCOE,NOCC,TCOE)
          !WL = WL+NOCC*(VINT_CI(LIST+1)+TCOE*VINT_CI(LIST))
          WL = VLOP0*NOCC
          call TRANS_IJKL_INTPOS(LRA,LRI,LRK,LRK,NXO)
          call PRODAB_2(3,JPEL,IWDL1,IWDR,JWL,JWR,WL,JPER,NXO)
          WL = VLOP0*NOCC*TCOE
          call TRANS_IJKL_INTPOS(LRA,LRK,LRI,LRK,NXO)
          call PRODAB_2(3,JPEL,IWDL1,IWDR,JWL,JWR,WL,JPER,NXO)

        end do
        !WL = WL*VLOP0
        ! SD(6-15) D&r&l(33)B^l(13)C'(22)
        VLOP0 = W0*W0SD15
        do LRK=1,LRI-1
          !LIST = LIST3(LRI,LRA,LRK)
          !WL = WL-VLOP0*(2*VINT_CI(LIST+1)-vint_ci(LIST))
          WL = -VLOP0*2
          call TRANS_IJKL_INTPOS(LRA,LRI,LRK,LRK,NXO)
          call PRODAB_2(3,JPEL,IWDL1,IWDR,JWL,JWR,WL,JPER,NXO)
          WL = VLOP0
          call TRANS_IJKL_INTPOS(LRA,LRK,LRI,LRK,NXO)
          call PRODAB_2(3,JPEL,IWDL1,IWDR,JWL,JWR,WL,JPER,NXO)

        end do
        !call PRODAB(3,JPEL,IWDL1,IWDR,JWL,JWR,WL,JPER)
      end if
      if ((LMIJ == JML) .and. (LMI == JMR)) then
        IWDR = JUD(LRI)
        IWDL1 = JUST(LRJ,LRI)
        W1SD10 = W1_SD(10)
        if (mod(NORB_DZ-LRJ,2) == 1) then
          W1SD10 = -W1SD10
        end if
        ! SD(6-10) D&r&l(12)B^l(23)
        VLOP1 = W1*W1SD10
        !LIST = LIST3(LRJ,LRA,LRI)
        !WL = -VLOP1*vint_ci(LIST)      !4.3
        WL = -VLOP1
        call TRANS_IJKL_INTPOS(LRA,LRI,LRJ,LRI,NXO)
        call PRODAB_2(3,JPEL,IWDL1,IWDR,JWL,JWR,WL,JPER,NXO)
      end if
    end if
    ! SD(6-4) A&r(23)-C'(12)-
    if (LMJ == JMR) then
      IWDR = JUD(LRJ)
      W0SD4 = W0_SD(4)
      W0SD16 = W0_SD(16)
      NI = mod(NORB_DZ-LRI,2)
      if (NI == 0) W0SD4 = -W0SD4
      if (NI == 0) W0SD16 = -W0SD16
      VLOP0 = W0*W0SD4
      !LIST = LIST3(LRI,LRA,LRI)
      !WL = VOINT(LRI,LRA)+VINT_CI(LIST)             !310,act_coe,61
      WL = VLOP0
      call PRODAB_1(3,JPEL,IWDL,IWDR,JWL,JWR,WL,JPER,LRI,LRA)
      call TRANS_IJKL_INTPOS(LRA,LRI,LRI,LRI,NXO)
      call PRODAB_2(3,JPEL,IWDL,IWDR,JWL,JWR,WL,JPER,NXO)

      !LIST = LIST3(LRI,LRA,LRJ)
      !WL = WL+VINT_CI(LIST+1)-(JB_SYS+2)*VINT_CI(LIST)
      WL = VLOP0
      call TRANS_IJKL_INTPOS(LRA,LRI,LRJ,LRJ,NXO)
      call PRODAB_2(3,JPEL,IWDL,IWDR,JWL,JWR,WL,JPER,NXO)
      WL = -VLOP0*(JB_SYS+2)
      call TRANS_IJKL_INTPOS(LRA,LRJ,LRI,LRJ,NXO)
      call PRODAB_2(3,JPEL,IWDL,IWDR,JWL,JWR,WL,JPER,NXO)

      do LR=LRI+1,NORB_DZ
        if (LR == LRJ) cycle
        !LIST = LIST3(LRI,LRA,LR)
        !WL = WL+2*VINT_CI(LIST+1)-VINT_CI(LIST)       !310:NEOC=2,C
        WL = VLOP0*2
        call TRANS_IJKL_INTPOS(LRA,LRI,LR,LR,NXO)
        call PRODAB_2(3,JPEL,IWDL,IWDR,JWL,JWR,WL,JPER,NXO)
        WL = -VLOP0
        call TRANS_IJKL_INTPOS(LRA,LR,LRI,LR,NXO)
        call PRODAB_2(3,JPEL,IWDL,IWDR,JWL,JWR,WL,JPER,NXO)

      end do
      do LRK=norb_dz+1,LRA
        !LIST = LIST3(LRI,LRA,LRK)
        KCOE = LPCOE(LRK)
        call NEOC(KCOE,NOCC,TCOE)
        !WL = WL+NOCC*(VINT_CI(LIST+1)+TCOE*VINT_CI(LIST))
        WL = VLOP0*NOCC
        call TRANS_IJKL_INTPOS(LRA,LRI,LRK,LRK,NXO)
        call PRODAB_2(3,JPEL,IWDL,IWDR,JWL,JWR,WL,JPER,NXO)
        WL = VLOP0*NOCC*TCOE
        call TRANS_IJKL_INTPOS(LRA,LRK,LRI,LRK,NXO)
        call PRODAB_2(3,JPEL,IWDL,IWDR,JWL,JWR,WL,JPER,NXO)

      end do
      !WL = WL*VLOP0
      ! SD(6-16) D&r&l(33)B^l(23)C'(12)
      VLOP0 = W0*W0SD16
      do LRK=1,LRI-1
        !LIST = LIST3(LRI,LRA,LRK)
        !WL = WL-VLOP0*(2*VINT_CI(LIST+1)-vint_ci(LIST))
        WL = -VLOP0*2
        call TRANS_IJKL_INTPOS(LRA,LRI,LRK,LRK,NXO)
        call PRODAB_2(3,JPEL,IWDL,IWDR,JWL,JWR,WL,JPER,NXO)
        WL = VLOP0
        call TRANS_IJKL_INTPOS(LRA,LRK,LRI,LRK,NXO)
        call PRODAB_2(3,JPEL,IWDL,IWDR,JWL,JWR,WL,JPER,NXO)

      end do
      !call PRODAB(3,JPEL,IWDL,IWDR,JWL,JWR,WL,JPER)
    end if

    ! SD(6-5) A&r(23)B&r(13)B^r(32)
    do LRD=LRJ+1,NORB_DZ
      LMD = LSM_INN(LRD)
      if (LMD /= JMR) cycle
      IWDR = JUD(LRD)
      W0SD5 = W0_SD(5)
      W1SD5 = W1_SD(5)
      NI = mod(LRJ-LRI+NORB_DZ-LRD,2)
      if (NI == 0) W0SD5 = -W0SD5
      if (NI == 0) W1SD5 = -W1SD5
      VLOP0 = W0*W0SD5
      VLOP1 = W1*W1SD5
      !LIST = LIST4(LRI,LRJ,LRD,LRA)
      !WL = (VLOP0-VLOP1)*vint_ci(LIST+2)+(VLOP0+VLOP1)*vint_ci(LIST)  !1.3
      WL = VLOP0-VLOP1
      call TRANS_IJKL_INTPOS(LRA,LRI,LRJ,LRD,NXO)
      call PRODAB_2(3,JPEL,IWDL,IWDR,JWL,JWR,WL,JPER,NXO)
      WL = VLOP0+VLOP1
      call TRANS_IJKL_INTPOS(LRA,LRJ,LRD,LRI,NXO)
      call PRODAB_2(3,JPEL,IWDL,IWDR,JWL,JWR,WL,JPER,NXO)

      !call PRODAB(3,JPEL,IWDL,IWDR,JWL,JWR,WL,JPER)
    end do
    if (JB_SYS > 0) then
      ! SD(6-6) A&r(13)B&r(23)B^r(32)
      do LRD=LRJ+1,NORB_DZ
        LMD = LSM_INN(LRD)
        if (LMD /= JMR) cycle
        IWDR = JUD(LRD)
        W0SD6 = W0_SD(6)
        W1SD6 = W1_SD(6)
        NI = mod(LRJ-LRI+NORB_DZ-LRD,2)
        if (NI == 0) W0SD6 = -W0SD6
        if (NI == 0) W1SD6 = -W1SD6
        VLOP0 = W0*W0SD6
        VLOP1 = W1*W1SD6
        !LIST = LIST4(LRI,LRJ,LRD,LRA)
        !WL = (VLOP0-VLOP1)*vint_ci(LIST+2)+(VLOP0+VLOP1)*vint_ci(LIST)   !1.3
        WL = VLOP0-VLOP1
        call TRANS_IJKL_INTPOS(LRA,LRI,LRJ,LRD,NXO)
        call PRODAB_2(3,JPEL,IWDL1,IWDR,JWL,JWR,WL,JPER,NXO)
        WL = VLOP0+VLOP1
        call TRANS_IJKL_INTPOS(LRA,LRJ,LRD,LRI,NXO)
        call PRODAB_2(3,JPEL,IWDL1,IWDR,JWL,JWR,WL,JPER,NXO)

        !call PRODAB(3,JPEL,IWDL1,IWDR,JWL,JWR,WL,JPER)
      end do
      ! SD(6-7) A&r(13)B&l(32)B^l(23)
      do LRD=LRI+1,LRJ-1
        LMD = LSM_INN(LRD)
        if (LMD /= JMR) cycle
        IWDR = JUD(LRD)
        W0SD7 = W0_SD(7)
        W1SD7 = W1_SD(7)
        NI = mod(LRD-LRI+NORB_DZ-LRJ,2)
        if (NI == 0) W0SD7 = -W0SD7
        if (NI == 0) W1SD7 = -W1SD7
        VLOP0 = W0*W0SD7
        VLOP1 = W1*W1SD7
        !LIST = LIST4(LRI,LRD,LRJ,LRA)
        !WL = (VLOP0-VLOP1)*vint_ci(LIST+2)-2*VLOP0*vint_ci(LIST+1)   !1.2
        WL = VLOP0-VLOP1
        call TRANS_IJKL_INTPOS(LRA,LRI,LRD,LRJ,NXO)
        call PRODAB_2(3,JPEL,IWDL1,IWDR,JWL,JWR,WL,JPER,NXO)
        WL = -2*VLOP0
        call TRANS_IJKL_INTPOS(LRA,LRJ,LRD,LRI,NXO)
        call PRODAB_2(3,JPEL,IWDL1,IWDR,JWL,JWR,WL,JPER,NXO)
        !call PRODAB(3,JPEL,IWDL1,IWDR,JWL,JWR,WL,JPER)
      end do
    end if
    ! SD(6-8) A&r(23)B&l(32)B^l(13)
    do LRD=LRI+1,LRJ-1
      LMD = LSM_INN(LRD)
      if (LMD /= JMR) cycle
      IWDR = JUD(LRD)
      W0SD8 = W0_SD(8)
      W1SD8 = W1_SD(8)
      NI = mod(LRD-LRI+NORB_DZ-LRJ,2)
      if (NI == 0) W0SD8 = -W0SD8
      if (NI == 0) W1SD8 = -W1SD8
      VLOP0 = W0*W0SD8
      VLOP1 = W1*W1SD8
      !LIST = LIST4(LRI,LRD,LRJ,LRA)
      !WL = (VLOP0-VLOP1)*vint_ci(LIST+2)-2*VLOP0*vint_ci(LIST+1)   !1.2
      WL = VLOP0-VLOP1
      call TRANS_IJKL_INTPOS(LRA,LRI,LRD,LRJ,NXO)
      call PRODAB_2(3,JPEL,IWDL,IWDR,JWL,JWR,WL,JPER,NXO)
      WL = -2*VLOP0
      call TRANS_IJKL_INTPOS(LRA,LRJ,LRD,LRI,NXO)
      call PRODAB_2(3,JPEL,IWDL,IWDR,JWL,JWR,WL,JPER,NXO)

      !call PRODAB(3,JPEL,IWDL,IWDR,JWL,JWR,WL,JPER)
    end do
  end do
end do

return

end subroutine SD_HEAD_DBL_TAIL_ACT_G

subroutine SDD_HEAD_DBL_TAIL_ACT_G(LRA,LPCOE)
!**********************************************
!     LRA, ....partial loop.......
!**********************************************

use gugaci_global, only: jb_sys, jml, jmr, jpel, jper, jud, just, jwl, jwr, lsm_inn, norb_dz, norb_frz, norb_inn, w0, w0_sd1, w1, &
                         w1_sd1
use Symmetry_Info, only: Mul
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: lra, lpcoe(norb_dz+1:norb_inn)
integer(kind=iwp) :: iwdl, iwdl1, iwdr, jmlr, kcoe, lmd, lmi, lmij, lmj, lr, lrd, lri, lrj, lrk, ni, nocc, nxo
real(kind=wp) :: tcoe, vlop0, vlop1, w0sd1, w0sd10, w0sd11, w0sd12, w0sd13, w0sd2, w0sd3, w0sd4, w0sd5, w0sd6, w0sd7, w0sd8, &
                 w0sd9, w1sd10, w1sd2, w1sd3, w1sd4, w1sd5, w1sd6, w1sd7, w1sd8, wl

JMLR = Mul(JML,JMR)
! SD1(8-1) A&r(01)-
! SD1(8-2) C(11)A&r(23)-
! SD1(8-3) A&r(13)C'(21)-
! SD1(8-4) A&r(23)C'(11)-
! SD1(8-5) A&r(13)B&r(23)B^r(31)
! SD1(8-6) A&r(23)B&r(13)B^r(31)
! SD1(8-7) A&r(13)B&l(31)B^l(23)
! SD1(8-8) A&r(23)B&l(31)B^l(13)
! SD1(8-9) D&r&r(03)B^r(31)
! SD1(8-10) D&r&l(11)B^l(23)
! SD1(8-11) D&r&l(33)B^l(01)
! SD1(8-12) D&r&l(33)B^l(23)
! SD1(8-13) D&r&l(33)C"(13)B^l(23)
! SD1(8-14) D&r&l(33)B^l(11)C'(23)
! SD1(8-15) D&r&l(33)B^l(23)C'(11)
if (JML == 1) then
  do LRI=NORB_FRZ+1,NORB_DZ
    ! SD1(8-1) A&r(01)-
    LMI = LSM_INN(LRI)
    IWDL = JUST(LRI,LRI)
    W0SD1 = W0_SD1(1)
    W0SD11 = W0_SD1(9)
    NI = mod(NORB_DZ-LRI,2)
    if (NI == 1) then
      W0SD1 = -W0SD1
      W0SD11 = -W0SD11
    end if
    if (LMI == JMLR) then
      IWDR = JUD(LRI)
      VLOP0 = W0*W0SD1
      !WL = VLOP0*VOINT(LRI,LRA)
      WL = VLOP0
      call PRODAB_1(3,JPEL,IWDL,IWDR,JWL,JWR,WL,JPER,LRI,LRA)

      do LR=LRI+1,NORB_DZ
        !LIST = LIST3(LRI,LRA,LR)
        !WL = WL+VLOP0*(2*VINT_CI(LIST+1)-VINT_CI(LIST)) !  310:NEOC=2
        WL = VLOP0*2
        call TRANS_IJKL_INTPOS(LRA,LRI,LR,LR,NXO)
        call PRODAB_2(3,JPEL,IWDL,IWDR,JWL,JWR,WL,JPER,NXO)
        WL = -VLOP0
        call TRANS_IJKL_INTPOS(LRA,LR,LRI,LR,NXO)
        call PRODAB_2(3,JPEL,IWDL,IWDR,JWL,JWR,WL,JPER,NXO)

      end do
      do LRK=norb_dz+1,LRA
        !LIST = LIST3(LRI,LRA,LRK)
        KCOE = LPCOE(LRK)
        call NEOC(KCOE,NOCC,TCOE)
        !WL = WL+VLOP0*NOCC*(VINT_CI(LIST+1)+TCOE*VINT_CI(LIST))
        WL = VLOP0*NOCC
        call TRANS_IJKL_INTPOS(LRA,LRI,LRK,LRK,NXO)
        call PRODAB_2(3,JPEL,IWDL,IWDR,JWL,JWR,WL,JPER,NXO)
        WL = VLOP0*NOCC*TCOE
        call TRANS_IJKL_INTPOS(LRA,LRK,LRI,LRK,NXO)
        call PRODAB_2(3,JPEL,IWDL,IWDR,JWL,JWR,WL,JPER,NXO)

      end do
      !WL = WL*VLOP0
      ! SD1(8-11) D&rl(33)B^l(01)
      VLOP0 = W0*W0SD11
      do LRK=1,LRI-1
        !LIST = LIST3(LRI,LRA,LRK)
        !WL = WL+VLOP0*(vint_ci(LIST)-2*VINT_CI(LIST+1))
        WL = VLOP0
        call TRANS_IJKL_INTPOS(LRA,LRK,LRI,LRK,NXO)
        call PRODAB_2(3,JPEL,IWDL,IWDR,JWL,JWR,WL,JPER,NXO)
        WL = -VLOP0*2
        call TRANS_IJKL_INTPOS(LRA,LRI,LRK,LRK,NXO)
        call PRODAB_2(3,JPEL,IWDL,IWDR,JWL,JWR,WL,JPER,NXO)

      end do
      !call PRODAB(3,JPEL,IWDL,IWDR,JWL,JWR,WL,JPER)
    end if

    ! SD1(8-9) D&r&r(03)B^r(31)
    do LRD=LRI+1,NORB_DZ
      LMD = LSM_INN(LRD)
      if (LMD /= JMR) cycle
      W0SD9 = W0_SD1(9)
      NI = mod(NORB_DZ-LRD,2)
      if (NI == 1) W0SD9 = -W0SD9
      IWDR = JUD(LRD)
      VLOP0 = W0*W0SD9
      !LIST = LIST3(LRD,LRA,LRI)
      !WL = VINT_CI(LIST)*VLOP0
      WL = VLOP0
      call TRANS_IJKL_INTPOS(LRA,LRI,LRD,LRI,NXO)
      call PRODAB_2(3,JPEL,IWDL,IWDR,JWL,JWR,WL,JPER,NXO)

      !call PRODAB(3,JPEL,IWDL,IWDR,JWL,JWR,WL,JPER)
    end do
  end do
end if

do LRI=NORB_FRZ+1,NORB_DZ
  LMI = LSM_INN(LRI)
  do LRJ=LRI+1,NORB_DZ
    if (LRI == LRJ) cycle
    LMJ = LSM_INN(LRJ)
    LMIJ = Mul(LMI,LMJ)
    if (LMIJ /= JML) cycle
    W0SD2 = W0_SD1(2)
    W1SD2 = W1_SD1(2)
    W0SD3 = W0_SD1(3)
    W1SD3 = W1_SD1(3)
    W0SD4 = W0_SD1(4)
    W1SD4 = W1_SD1(4)
    W0SD10 = W0_SD1(10)
    W1SD10 = W1_SD1(10)
    W0SD11 = W0_SD1(11)
    W0SD12 = W0_SD1(12)
    W0SD13 = W0_SD1(13)
    NI = mod(NORB_DZ-LRJ,2)
    if (NI == 1) then
      W0SD2 = -W0SD2
      W1SD2 = -W1SD2
      W0SD10 = -W0SD10
      W1SD10 = -W1SD10
      W0SD11 = -W0SD11
    end if
    if (mod(NORB_DZ-LRI,2) == 1) then
      W0SD3 = -W0SD3
      W1SD3 = -W1SD3
      W0SD4 = -W0SD4
      W1SD4 = -W1SD4
      W0SD12 = -W0SD12
      W0SD13 = -W0SD13
    end if
    IWDL = JUST(LRJ,LRI)
    IWDL1 = JUST(LRI,LRJ)
    !*******************************************************************
    ! SD1(8-2) C(11)-A&r(23)-
    if (LMI == JMR) then
      IWDL = JUST(LRJ,LRI)
      IWDR = JUD(LRI)
      VLOP0 = W0*W0SD2
      !LIST = LIST3(LRJ,LRA,LRJ)
      !WL = VLOP0*(VOINT(LRJ,LRA)+VINT_CI(LIST))
      WL = VLOP0
      call PRODAB_1(3,JPEL,IWDL,IWDR,JWL,JWR,WL,JPER,LRJ,LRA)
      call TRANS_IJKL_INTPOS(LRA,LRJ,LRJ,LRJ,NXO)
      call PRODAB_2(3,JPEL,IWDL,IWDR,JWL,JWR,WL,JPER,NXO)

      do LR=LRJ+1,NORB_DZ
        !LIST = LIST3(LRJ,LRA,LR)
        !WL = WL+VLOP0*(2*VINT_CI(LIST+1)-VINT_CI(LIST)) !  310:NEOC
        WL = VLOP0*2
        call TRANS_IJKL_INTPOS(LRA,LRJ,LR,LR,NXO)
        call PRODAB_2(3,JPEL,IWDL,IWDR,JWL,JWR,WL,JPER,NXO)
        WL = -VLOP0
        call TRANS_IJKL_INTPOS(LRA,LR,LRJ,LR,NXO)
        call PRODAB_2(3,JPEL,IWDL,IWDR,JWL,JWR,WL,JPER,NXO)

      end do
      do LRK=norb_dz+1,LRA
        !LIST = LIST3(LRJ,LRA,LRK)
        KCOE = LPCOE(LRK)
        call NEOC(KCOE,NOCC,TCOE)
        !WL = WL+VLOP0*NOCC*(VINT_CI(LIST+1)+TCOE*VINT_CI(LIST))
        WL = VLOP0*NOCC
        call TRANS_IJKL_INTPOS(LRA,LRJ,LRK,LRK,NXO)
        call PRODAB_2(3,JPEL,IWDL,IWDR,JWL,JWR,WL,JPER,NXO)
        WL = VLOP0*NOCC*TCOE
        call TRANS_IJKL_INTPOS(LRA,LRK,LRJ,LRK,NXO)
        call PRODAB_2(3,JPEL,IWDL,IWDR,JWL,JWR,WL,JPER,NXO)

      end do
      !WL = WL*VLOP0
      ! SD1(8-10) D&r&l(11)B^l(23)
      VLOP0 = W0*W0SD10
      VLOP1 = W1*W1SD10
      !LIST = LIST3(LRJ,LRA,LRI)
      !WL = WL+(VLOP0-VLOP1)*vint_ci(LIST)-2*VLOP0*VINT_CI(LIST+1)

      WL = VLOP0-VLOP1
      call TRANS_IJKL_INTPOS(LRA,LRI,LRJ,LRI,NXO)
      call PRODAB_2(3,JPEL,IWDL,IWDR,JWL,JWR,WL,JPER,NXO)
      WL = -VLOP0*2
      call TRANS_IJKL_INTPOS(LRA,LRJ,LRI,LRI,NXO)
      call PRODAB_2(3,JPEL,IWDL,IWDR,JWL,JWR,WL,JPER,NXO)

      ! SD1(8-11) (11)D&r&l(33)B^l(23)
      ! SD1(8-11) D&r&l(33)C"(11)B^l(23)
      VLOP0 = W0*W0SD11
      do LRK=1,LRJ-1
        if (LRK == LRI) cycle
        !LIST = LIST3(LRJ,LRA,LRK)
        !WL = WL+VLOP0*(vint_ci(LIST)-2*VINT_CI(LIST+1))
        WL = VLOP0
        call TRANS_IJKL_INTPOS(LRA,LRK,LRJ,LRK,NXO)
        call PRODAB_2(3,JPEL,IWDL,IWDR,JWL,JWR,WL,JPER,NXO)
        WL = -VLOP0*2
        call TRANS_IJKL_INTPOS(LRA,LRJ,LRK,LRK,NXO)
        call PRODAB_2(3,JPEL,IWDL,IWDR,JWL,JWR,WL,JPER,NXO)

      end do
      !call PRODAB(3,JPEL,IWDL,IWDR,JWL,JWR,WL,JPER)
    end if
    ! SD1(8-3) A&r(13)-C'(21)-
    if (LMJ == JMR) then
      VLOP0 = -W0*W0SD3
      !-----------------------------------------------------------------
      ! lyb
      IWDL = JUST(LRJ,LRI)
      IWDR = JUD(LRJ)

      !LIST = LIST3(LRI,LRA,LRI)
      !WL = VOINT(LRI,LRA)+VINT_CI(LIST)
      WL = VLOP0
      !write(nf2,'(a3,5I8)') 'sdd',JPEL,IWDL,IWDR,JWL,JWR
      call PRODAB_1(3,JPEL,IWDL,IWDR,JWL,JWR,WL,JPER,LRI,LRA)
      call TRANS_IJKL_INTPOS(LRA,LRI,LRI,LRI,NXO)
      call PRODAB_2(3,JPEL,IWDL,IWDR,JWL,JWR,WL,JPER,NXO)

      !LIST = LIST3(LRI,LRA,LRJ)
      !WL = WL+VINT_CI(LIST+1)+JB_SYS*VINT_CI(LIST)
      WL = VLOP0
      call TRANS_IJKL_INTPOS(LRA,LRI,LRJ,LRJ,NXO)
      call PRODAB_2(3,JPEL,IWDL,IWDR,JWL,JWR,WL,JPER,NXO)
      WL = VLOP0*JB_SYS
      call TRANS_IJKL_INTPOS(LRA,LRJ,LRI,LRJ,NXO)
      call PRODAB_2(3,JPEL,IWDL,IWDR,JWL,JWR,WL,JPER,NXO)

      do LR=LRI+1,NORB_DZ
        if (LR == LRJ) cycle
        !LIST = LIST3(LRI,LRA,LR)
        !WL = WL+2*VINT_CI(LIST+1)-VINT_CI(LIST)       !310:NEOC=2,C
        WL = VLOP0*2
        call TRANS_IJKL_INTPOS(LRA,LRI,LR,LR,NXO)
        call PRODAB_2(3,JPEL,IWDL,IWDR,JWL,JWR,WL,JPER,NXO)
        WL = -VLOP0
        call TRANS_IJKL_INTPOS(LRA,LR,LRI,LR,NXO)
        call PRODAB_2(3,JPEL,IWDL,IWDR,JWL,JWR,WL,JPER,NXO)

      end do
      do LRK=norb_dz+1,LRA
        !LIST = LIST3(LRI,LRA,LRK)
        KCOE = LPCOE(LRK)
        call NEOC(KCOE,NOCC,TCOE)
        !WL = WL+NOCC*(VINT_CI(LIST+1)+TCOE*VINT_CI(LIST))
        WL = VLOP0*NOCC
        call TRANS_IJKL_INTPOS(LRA,LRI,LRK,LRK,NXO)
        call PRODAB_2(3,JPEL,IWDL,IWDR,JWL,JWR,WL,JPER,NXO)
        WL = VLOP0*NOCC*TCOE
        call TRANS_IJKL_INTPOS(LRA,LRK,LRI,LRK,NXO)
        call PRODAB_2(3,JPEL,IWDL,IWDR,JWL,JWR,WL,JPER,NXO)
      end do
      !WL = WL*VLOP0
      ! SD1(8-12) Drl(33)-BL(13)-C'(21)-
      !IWDL = JUST(LRJ,LRI)
      !IWDR = JUD(LRJ)
      VLOP0 = -W0*W0SD12
      do LRK=1,LRI-1
        !LIST = LIST3(LRI,LRA,LRK)
        !WL = WL-VLOP0*(2*VINT_CI(LIST+1)-vint_ci(LIST))
        WL = -VLOP0*2
        call TRANS_IJKL_INTPOS(LRA,LRI,LRK,LRK,NXO)
        call PRODAB_2(3,JPEL,IWDL,IWDR,JWL,JWR,WL,JPER,NXO)
        WL = VLOP0
        call TRANS_IJKL_INTPOS(LRA,LRK,LRI,LRK,NXO)
        call PRODAB_2(3,JPEL,IWDL,IWDR,JWL,JWR,WL,JPER,NXO)

      end do
      !call PRODAB(3,JPEL,IWDL,IWDR,JWL,JWR,WL,JPER)
    end if
    ! SD1(8-4) A&r(23)-C'(11)-
    if (LMJ == JMR) then
      !-----------------------------------------------------------------
      ! lyb
      IWDL = JUST(LRI,LRJ)
      IWDR = JUD(LRJ)
      VLOP0 = -W0*W0SD4
      !LIST = LIST3(LRI,LRA,LRI)
      !WL = VOINT(LRI,LRA)+VINT_CI(LIST)
      WL = VLOP0
      call PRODAB_1(3,JPEL,IWDL,IWDR,JWL,JWR,WL,JPER,LRI,LRA)
      call TRANS_IJKL_INTPOS(LRA,LRI,LRI,LRI,NXO)
      call PRODAB_2(3,JPEL,IWDL,IWDR,JWL,JWR,WL,JPER,NXO)

      !LIST = LIST3(LRI,LRA,LRJ)
      !WL = WL+VINT_CI(LIST+1)-VINT_CI(LIST)
      WL = VLOP0
      call TRANS_IJKL_INTPOS(LRA,LRI,LRJ,LRJ,NXO)
      call PRODAB_2(3,JPEL,IWDL,IWDR,JWL,JWR,WL,JPER,NXO)
      WL = -VLOP0
      call TRANS_IJKL_INTPOS(LRA,LRJ,LRI,LRJ,NXO)
      call PRODAB_2(3,JPEL,IWDL,IWDR,JWL,JWR,WL,JPER,NXO)

      do LR=LRI+1,NORB_DZ
        if (LR == LRJ) cycle
        !LIST = LIST3(LRI,LRA,LR)
        !WL = WL+2*VINT_CI(LIST+1)-VINT_CI(LIST)       !310:NEOC=2,C
        WL = VLOP0*2
        call TRANS_IJKL_INTPOS(LRA,LRI,LR,LR,NXO)
        call PRODAB_2(3,JPEL,IWDL,IWDR,JWL,JWR,WL,JPER,NXO)
        WL = -VLOP0
        call TRANS_IJKL_INTPOS(LRA,LR,LRI,LR,NXO)
        call PRODAB_2(3,JPEL,IWDL,IWDR,JWL,JWR,WL,JPER,NXO)

      end do
      do LRK=norb_dz+1,LRA
        !LIST = LIST3(LRI,LRA,LRK)
        KCOE = LPCOE(LRK)
        call NEOC(KCOE,NOCC,TCOE)
        !WL = WL+NOCC*(VINT_CI(LIST+1)+TCOE*VINT_CI(LIST))
        WL = VLOP0*NOCC
        call TRANS_IJKL_INTPOS(LRA,LRI,LRK,LRK,NXO)
        call PRODAB_2(3,JPEL,IWDL,IWDR,JWL,JWR,WL,JPER,NXO)
        WL = VLOP0*NOCC*TCOE
        call TRANS_IJKL_INTPOS(LRA,LRK,LRI,LRK,NXO)
        call PRODAB_2(3,JPEL,IWDL,IWDR,JWL,JWR,WL,JPER,NXO)

      end do
      !WL = WL*VLOP0
      ! SD1(8-13) Drl(33)-BL(23)-C'(11)-
      !IWDL = JUST(LRI,LRJ)
      !IWDR = JUD(LRJ)
      VLOP0 = -W0*W0SD13
      do LRK=1,LRI-1
        !LIST = LIST3(LRI,LRA,LRK)
        !WL = WL-VLOP0*(2*VINT_CI(LIST+1)-vint_ci(LIST))
        WL = -VLOP0*2
        call TRANS_IJKL_INTPOS(LRA,LRI,LRK,LRK,NXO)
        call PRODAB_2(3,JPEL,IWDL,IWDR,JWL,JWR,WL,JPER,NXO)
        WL = VLOP0
        call TRANS_IJKL_INTPOS(LRA,LRK,LRI,LRK,NXO)
        call PRODAB_2(3,JPEL,IWDL,IWDR,JWL,JWR,WL,JPER,NXO)

      end do
      !call PRODAB(3,JPEL,IWDL,IWDR,JWL,JWR,WL,JPER)
    end if
    ! SD1(8-5) A&r(13)B&r(23)B^r(31)
    do LRD=LRJ+1,NORB_DZ
      LMD = LSM_INN(LRD)
      if (LMD /= JMR) cycle
      IWDL = JUST(LRJ,LRI)
      IWDR = JUD(LRD)
      W0SD5 = W0_SD1(5)
      W1SD5 = W1_SD1(5)
      NI = mod(LRJ-LRI+NORB_DZ-LRD,2)
      if (NI == 0) W0SD5 = -W0SD5
      if (NI == 0) W1SD5 = -W1SD5
      VLOP0 = W0*W0SD5
      VLOP1 = W1*W1SD5
      !LIST = LIST4(LRI,LRJ,LRD,LRA)
      !WL = (VLOP0-VLOP1)*vint_ci(LIST+2)+(VLOP0+VLOP1)*vint_ci(LIST)      !1.3
      WL = VLOP0-VLOP1
      call TRANS_IJKL_INTPOS(LRA,LRI,LRJ,LRD,NXO)
      call PRODAB_2(3,JPEL,IWDL,IWDR,JWL,JWR,WL,JPER,NXO)
      WL = VLOP0+VLOP1
      call TRANS_IJKL_INTPOS(LRA,LRJ,LRD,LRI,NXO)
      call PRODAB_2(3,JPEL,IWDL,IWDR,JWL,JWR,WL,JPER,NXO)
      !call PRODAB(3,JPEL,IWDL,IWDR,JWL,JWR,WL,JPER)
      ! SD1(8-6) A&r(23)B&r(13)B^r(31)
      IWDL = JUST(LRI,LRJ)
      W0SD6 = W0_SD1(6)
      W1SD6 = W1_SD1(6)
      NI = mod(LRJ-LRI+NORB_DZ-LRD,2)
      if (NI == 0) W0SD6 = -W0SD6
      if (NI == 0) W1SD6 = -W1SD6
      VLOP0 = W0*W0SD6
      VLOP1 = W1*W1SD6
      !LIST = LIST4(LRI,LRJ,LRD,LRA)
      !WL = (VLOP0-VLOP1)*vint_ci(LIST+2)+(VLOP0+VLOP1)*vint_ci(LIST)   !1.3
      WL = VLOP0-VLOP1
      call TRANS_IJKL_INTPOS(LRA,LRI,LRJ,LRD,NXO)
      call PRODAB_2(3,JPEL,IWDL1,IWDR,JWL,JWR,WL,JPER,NXO)
      WL = VLOP0+VLOP1
      call TRANS_IJKL_INTPOS(LRA,LRJ,LRD,LRI,NXO)
      call PRODAB_2(3,JPEL,IWDL1,IWDR,JWL,JWR,WL,JPER,NXO)

      !call PRODAB(3,JPEL,IWDL1,IWDR,JWL,JWR,WL,JPER)
    end do
    ! SD1(8-7) A&r(13)B&l(31)B^l(23)
    do LRD=LRI+1,LRJ-1
      LMD = LSM_INN(LRD)
      if (LMD /= JMR) cycle
      IWDR = JUD(LRD)
      IWDL = JUST(LRJ,LRI)
      W0SD7 = W0_SD1(7)
      W1SD7 = W1_SD1(7)
      NI = mod(LRD-LRI+NORB_DZ-LRJ,2)
      if (NI == 0) W0SD7 = -W0SD7
      if (NI == 0) W1SD7 = -W1SD7
      VLOP0 = W0*W0SD7
      VLOP1 = W1*W1SD7
      !LIST = LIST4(LRI,LRD,LRJ,LRA)
      !WL = (VLOP0-VLOP1)*vint_ci(LIST+2)-2*VLOP0*vint_ci(LIST+1)      !1.2
      WL = VLOP0-VLOP1
      call TRANS_IJKL_INTPOS(LRA,LRI,LRD,LRJ,NXO)
      call PRODAB_2(3,JPEL,IWDL,IWDR,JWL,JWR,WL,JPER,NXO)
      WL = -2*VLOP0
      call TRANS_IJKL_INTPOS(LRA,LRJ,LRD,LRI,NXO)
      call PRODAB_2(3,JPEL,IWDL,IWDR,JWL,JWR,WL,JPER,NXO)

      !call PRODAB(3,JPEL,IWDL,IWDR,JWL,JWR,WL,JPER)
      ! SD1(8-8) A&r(23)B&l(31)B^l(13)
      IWDL = JUST(LRI,LRJ)
      W0SD8 = W0_SD1(8)
      W1SD8 = W1_SD1(8)
      NI = mod(LRD-LRI+NORB_DZ-LRJ,2)
      if (NI == 0) W0SD8 = -W0SD8
      if (NI == 0) W1SD8 = -W1SD8
      VLOP0 = W0*W0SD8
      VLOP1 = W1*W1SD8
      !LIST = LIST4(LRI,LRD,LRJ,LRA)
      !WL = (VLOP0-VLOP1)*vint_ci(LIST+2)-2*VLOP0*vint_ci(LIST+1)   !1.2
      WL = VLOP0-VLOP1
      call TRANS_IJKL_INTPOS(LRA,LRI,LRD,LRJ,NXO)
      call PRODAB_2(3,JPEL,IWDL1,IWDR,JWL,JWR,WL,JPER,NXO)
      WL = -2*VLOP0
      call TRANS_IJKL_INTPOS(LRA,LRJ,LRD,LRI,NXO)
      call PRODAB_2(3,JPEL,IWDL1,IWDR,JWL,JWR,WL,JPER,NXO)

      !call PRODAB(3,JPEL,IWDL1,IWDR,JWL,JWR,WL,JPER)
    end do
  end do
end do

return

end subroutine SDD_HEAD_DBL_TAIL_ACT_G

!subroutine TD_HEAD_DBL_TAIL_ACT_G(LRA,LPCOE)
!!**********************************************
!!     LRA, ....partial loop.......
!!**********************************************
!
!use gugaci_global, only: jml, jmr, jpel, jper, jud, just, jwl, jwr, lsm_inn, norb_dz, norb_frz, norb_inn, w0, w0_td, w1, w1_td
!use Symmetry_Info, only: Mul
!use Definitions, only: wp, iwp
!
!implicit none
!integer(kind=iwp), intent(in) :: lra, lpcoe(norb_dz+1:norb_inn)
!integer(kind=iwp) :: iwdl, iwdr, jmlr, kcoe, lmd, lmi, lmij, lmj, lr, lrd, lri, lrj, lrk, ni, nocc, nxo
!real(kind=wp) :: tcoe, vlop0, vlop1, w0td1, w0td2, w0td3, w0td4, w0td5, w1td2, w1td3, wl
!
!JMLR = Mul(JML,JMR)
!! TD(13-1) (22)A&(23)
!! TD(13-1) A&(23)C'(22)
!! TD(13-5) (22)D&&l(33)B^l(23)
!do LRI=NORB_FRZ+1,NORB_DZ
!  LMI = LSM_INN(LRI)
!  if (LMI /= JMLR) cycle
!  W0TD1 = W0_TD(1)
!  W0TD4 = W0_TD(4)
!  W0TD5 = W0_TD(5)
!  NI = mod(NORB_DZ-LRI,2)
!  if (NI == 1) W0TD1 = -W0TD1
!  if (NI == 1) W0TD4 = -W0TD4
!  if (NI == 1) W0TD5 = -W0TD5
!
!  ! TD(13-1) A&(23)C'(22)
!  do LRD=LRI+1,NORB_DZ
!    LMD = LSM_INN(LRD)
!    if (LMD /= JMR) cycle
!    IWDL = JUST(LRI,LRD)
!    IWDR = JUD(LRD)
!    VLOP0 = -W0*W0TD1
!    !LIST = LIST3(LRI,LRA,LRI)
!    !WL = VOINT(LRI,LRA)+VINT_CI(LIST)             !310,act_coe,610,
!    WL = VLOP0
!    call PRODAB_1(3,JPEL,IWDL,IWDR,JWL,JWR,WL,JPER,LRI,LRA)
!    call TRANS_IJKL_INTPOS(LRA,LRI,LRI,LRI,NXO)
!    call PRODAB_2(3,JPEL,IWDL,IWDR,JWL,JWR,WL,JPER,NXO)
!
!    !LIST = LIST3(LRI,LRA,LRD)
!    !WL = WL+VINT_CI(LIST+1)                          !310 C'(22) COE
!    WL = VLOP0
!    call TRANS_IJKL_INTPOS(LRA,LRD,LRI,LRD,NXO)
!    call PRODAB_2(3,JPEL,IWDL,IWDR,JWL,JWR,WL,JPER,NXO)
!
!    do LR=LRI+1,NORB_DZ
!      if (LR == LRD) cycle
!      !LIST = LIST3(LRI,LRA,LR)
!      !WL = WL+2*VINT_CI(LIST+1)-VINT_CI(LIST)       !310:NEOC=2,COE
!      WL = 2*VLOP0
!      call TRANS_IJKL_INTPOS(LRA,LRI,LR,LR,NXO)
!      call PRODAB_2(3,JPEL,IWDL,IWDR,JWL,JWR,WL,JPER,NXO)
!      WL = -VLOP0
!      call TRANS_IJKL_INTPOS(LRA,LR,LRI,LR,NXO)
!      call PRODAB_2(3,JPEL,IWDL,IWDR,JWL,JWR,WL,JPER,NXO)
!
!    end do
!    do LRK=norb_dz+1,LRA
!      !LIST = LIST3(LRI,LRA,LRK)
!      KCOE = LPCOE(LRK)
!      call NEOC(KCOE,NOCC,TCOE)
!      !WL = WL+NOCC*(VINT_CI(LIST+1)+TCOE*VINT_CI(LIST))
!      call TRANS_IJKL_INTPOS(LRA,LRI,LRK,LRK,NXO)
!      WL = VLOP0*NOCC
!      call PRODAB_2(3,JPEL,IWDL,IWDR,JWL,JWR,WL,JPER,NXO)
!      WL = VLOP0*NOCC*TCOE
!      call TRANS_IJKL_INTPOS(LRA,LRK,LRI,LRK,NXO)
!      call PRODAB_2(3,JPEL,IWDL,IWDR,JWL,JWR,WL,JPER,NXO)
!
!    end do
!    !WL = WL*VLOP0
!    ! TD(13-5) D&rl(33)B^l(23)C'(22)          !CC (22)D&rl(33)B^l(23)???
!    VLOP0 = -W0*W0TD5
!    do LRK=1,LRI-1
!      !LIST = LIST3(LRI,LRA,LRK)
!      !WL = WL-VLOP0*(2*VINT_CI(LIST+1)-vint_ci(LIST))
!      WL = -2*VLOP0
!      call TRANS_IJKL_INTPOS(LRA,LRI,LRK,LRK,NXO)
!      call PRODAB_2(3,JPEL,IWDL,IWDR,JWL,JWR,WL,JPER,NXO)
!      WL = VLOP0
!      call TRANS_IJKL_INTPOS(LRA,LRK,LRI,LRK,NXO)
!      call PRODAB_2(3,JPEL,IWDL,IWDR,JWL,JWR,WL,JPER,NXO)
!
!    end do
!    !call PRODAB(3,JPEL,IWDL,IWDR,JWL,JWR,WL,JPER)
!  end do
!  !---------------------------------------------------------------------
!  do LRD=NORB_FRZ+1,LRI-1
!    LMD = LSM_INN(LRD)
!    if (LMD /= JMR) cycle
!    IWDL = JUST(LRD,LRI)
!    IWDR = JUD(LRD)
!    ! TD(13-1) (22)A&(23)
!    VLOP0 = W0*W0TD1
!    !LIST = LIST3(LRI,LRA,LRI)
!    !WL = VLOP0*(VOINT(LRI,LRA)+VINT_CI(LIST))             !310,act_
!    WL = VLOP0
!    call PRODAB_1(3,JPEL,IWDL,IWDR,JWL,JWR,WL,JPER,LRI,LRA)
!    call TRANS_IJKL_INTPOS(LRA,LRI,LRI,LRI,NXO)
!    call PRODAB_2(3,JPEL,IWDL,IWDR,JWL,JWR,WL,JPER,NXO)
!
!    do LR=LRI+1,NORB_DZ
!      !LIST = LIST3(LRI,LRA,LR)
!      !WL = WL+VLOP0*(2*VINT_CI(LIST+1)-VINT_CI(LIST)) !  310:NEOC=2
!      WL = 2*VLOP0
!      call TRANS_IJKL_INTPOS(LRA,LRI,LR,LR,NXO)
!      call PRODAB_2(3,JPEL,IWDL,IWDR,JWL,JWR,WL,JPER,NXO)
!      WL = -VLOP0
!      call TRANS_IJKL_INTPOS(LRA,LR,LRI,LR,NXO)
!      call PRODAB_2(3,JPEL,IWDL,IWDR,JWL,JWR,WL,JPER,NXO)
!
!    end do
!    do LRK=norb_dz+1,LRA
!      !LIST = LIST3(LRI,LRA,LRK)
!      KCOE = LPCOE(LRK)
!      call NEOC(KCOE,NOCC,TCOE)
!      !WL = WL+VLOP0*NOCC*(VINT_CI(LIST+1)+TCOE*VINT_CI(LIST))
!      WL = VLOP0*NOCC
!      call TRANS_IJKL_INTPOS(LRA,LRI,LRK,LRK,NXO)
!      call PRODAB_2(3,JPEL,IWDL,IWDR,JWL,JWR,WL,JPER,NXO)
!      WL = VLOP0*NOCC*TCOE
!      call TRANS_IJKL_INTPOS(LRA,LRK,LRI,LRK,NXO)
!      call PRODAB_2(3,JPEL,IWDL,IWDR,JWL,JWR,WL,JPER,NXO)
!    end do
!    !WL = WL*VLOP0
!    ! TD(13-4) D&r&l(22)B^l(23)
!    VLOP0 = W0*W0TD4
!    VLOP1 = W1*W0TD4
!    !LIST = LIST3(LRI,LRA,LRD)
!    !WL = WL+(VLOP0-VLOP1)*vint_ci(LIST)-2*VLOP0*VINT_CI(LIST+1)
!
!    WL = VLOP0-VLOP1
!    call TRANS_IJKL_INTPOS(LRA,LRD,LRI,LRD,NXO)
!    call PRODAB_2(3,JPEL,IWDL,IWDR,JWL,JWR,WL,JPER,NXO)
!    WL = -2*VLOP0
!    call TRANS_IJKL_INTPOS(LRA,LRI,LRD,LRD,NXO)
!    call PRODAB_2(3,JPEL,IWDL,IWDR,JWL,JWR,WL,JPER,NXO)
!    ! TD(13-5) D&rl(33)C"(22)B^l(23)
!    VLOP0 = W0*W0TD5
!    do LRK=1,LRI-1
!      if (LRK == LRD) cycle
!      !LIST = LIST3(LRI,LRA,LRK)
!      !WL = WL+VLOP0*(vint_ci(LIST)-2*VINT_CI(LIST+1))      !4.3
!      WL = VLOP0
!      call TRANS_IJKL_INTPOS(LRA,LRK,LRI,LRK,NXO)
!      call PRODAB_2(3,JPEL,IWDL,IWDR,JWL,JWR,WL,JPER,NXO)
!      WL = -2*VLOP0
!      call TRANS_IJKL_INTPOS(LRA,LRI,LRK,LRK,NXO)
!      call PRODAB_2(3,JPEL,IWDL,IWDR,JWL,JWR,WL,JPER,NXO)
!    end do
!    !call PRODAB(3,JPEL,IWDL,IWDR,JWL,JWR,WL,JPER)
!  end do
!end do
!do LRI=NORB_FRZ+1,NORB_DZ-1
!  LMI = LSM_INN(LRI)
!  do LRJ=LRI+1,NORB_DZ
!    LMJ = LSM_INN(LRJ)
!    LMIJ = Mul(LMI,LMJ)
!    if (LMIJ /= JML) cycle
!    IWDL = JUST(LRI,LRJ)
!
!    ! TD(13-2) A&(23)B&r(23)B^r(32)
!    do LRD=LRJ+1,NORB_DZ
!      LMD = LSM_INN(LRD)
!      if (LMD /= JMR) cycle
!      W0TD2 = W0_TD(2)
!      W1TD2 = W1_TD(2)
!      NI = mod(LRJ-LRI+NORB_DZ-LRD,2)
!      if (NI == 0) W0TD2 = -W0TD2
!      if (NI == 0) W1TD2 = -W1TD2
!
!      IWDR = JUD(LRD)
!      VLOP0 = W0*W0TD2
!      VLOP1 = W1*W1TD2
!      !LIST = LIST4(LRI,LRJ,LRD,LRA)
!      !WL = VLOP0*(vint_ci(LIST+2)+vint_ci(LIST))-VLOP1*(vint_ci(LIST+2)-vint_ci(LIST))  !1.3
!      WL = VLOP0-VLOP1
!      call TRANS_IJKL_INTPOS(LRA,LRI,LRJ,LRD,NXO)
!      call PRODAB_2(3,JPEL,IWDL,IWDR,JWL,JWR,WL,JPER,NXO)
!      WL = VLOP0+VLOP1
!      call TRANS_IJKL_INTPOS(LRA,LRJ,LRD,LRI,NXO)
!      call PRODAB_2(3,JPEL,IWDL,IWDR,JWL,JWR,WL,JPER,NXO)
!
!      !call PRODAB(3,JPEL,IWDL,IWDR,JWL,JWR,WL,JPER)
!    end do
!    ! TD(13-3) A&(23)B&l(32)B^l(23)
!    do LRD=LRI+1,LRJ-1
!      LMD = LSM_INN(LRD)
!      if (LMD /= JMR) cycle
!      IWDR = JUD(LRD)
!      W0TD3 = W0_TD(3)
!      W1TD3 = W1_TD(3)
!      NI = mod(LRD-LRI+NORB_DZ-LRJ,2)
!      if (NI == 0) W0TD3 = -W0TD3
!      if (NI == 0) W1TD3 = -W1TD3
!      VLOP0 = W0*W0TD3                !D6-8
!      VLOP1 = W1*W1TD3
!      !LIST = LIST4(LRI,LRD,LRJ,LRA)
!      !WL = VLOP0*(vint_ci(LIST+2)-2*vint_ci(LIST+1))-VLOP1*vint_ci(LIST+2)  !1.2
!      WL = VLOP0-VLOP1
!      call TRANS_IJKL_INTPOS(LRA,LRI,LRD,LRJ,NXO)
!      call PRODAB_2(3,JPEL,IWDL,IWDR,JWL,JWR,WL,JPER,NXO)
!      WL = -2*VLOP0
!      call TRANS_IJKL_INTPOS(LRA,LRJ,LRD,LRI,NXO)
!      call PRODAB_2(3,JPEL,IWDL,IWDR,JWL,JWR,WL,JPER,NXO)
!
!      !call PRODAB(3,JPEL,IWDL,IWDR,JWL,JWR,WL,JPER)
!    end do
!  end do
!end do
!
!return
!
!end subroutine TD_HEAD_DBL_TAIL_ACT_G

subroutine TTDD_HEAD_DBL_TAIL_ACT_G(LRA,LPCOE)
!**********************************************
!     LRA, ....partial loop.......
!**********************************************

use gugaci_global, only: jml, jmr, jpel, jper, jud, just, jwl, jwr, lsm_inn, norb_dz, norb_frz, norb_inn, w0, w0_t1d1, w1, w1_t1d1
use Symmetry_Info, only: Mul
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: lra, lpcoe(norb_dz+1:norb_inn)
integer(kind=iwp) :: iwdl, iwdr, jmlr, kcoe, lmd, lmi, lmij, lmj, lr, lrd, lri, lrj, lrk, ni, nocc, nxo
real(kind=wp) :: tcoe, vlop0, vlop1, w0td1, w0td2, w0td3, w0td4, w0td5, w1td2, w1td3, wl

JMLR = Mul(JML,JMR)
! T1D1(15-1) (11)A&(13)
! T1D1(15-1) A&(13)C'(11)
! T1D1(15-5) (11)D&&l(33)B^l(13)
do LRI=NORB_FRZ+1,NORB_DZ
  LMI = LSM_INN(LRI)
  if (LMI /= JMLR) cycle
  W0TD1 = W0_T1D1(1)
  W0TD4 = W0_T1D1(4)
  W0TD5 = W0_T1D1(5)
  NI = mod(NORB_DZ-LRI,2)
  if (NI == 1) then
    W0TD1 = -W0TD1
    W0TD4 = -W0TD4
    W0TD5 = -W0TD5
  end if
  ! T1D1(15-2) A&(13)C'(11)
  do LRD=LRI+1,NORB_DZ
    LMD = LSM_INN(LRD)
    if (LMD /= JMR) cycle
    IWDL = JUST(LRI,LRD)
    IWDR = JUD(LRD)
    VLOP0 = -W0*W0TD1
    !LIST = LIST3(LRI,LRA,LRI)
    !WL = VOINT(LRI,LRA)+VINT_CI(LIST)             !310,act_coe,610,
    WL = VLOP0
    call PRODAB_1(3,JPEL,IWDL,IWDR,JWL,JWR,WL,JPER,LRI,LRA)
    call TRANS_IJKL_INTPOS(LRA,LRI,LRI,LRI,NXO)
    call PRODAB_2(3,JPEL,IWDL,IWDR,JWL,JWR,WL,JPER,NXO)

    !LIST = LIST3(LRI,LRA,LRD)
    !WL = WL+VINT_CI(LIST+1)                          !310 C'(22) CO
    WL = VLOP0
    call TRANS_IJKL_INTPOS(LRA,LRI,LRD,LRD,NXO)
    call PRODAB_2(3,JPEL,IWDL,IWDR,JWL,JWR,WL,JPER,NXO)

    do LR=LRI+1,NORB_DZ
      if (LR == LRD) cycle
      !LIST = LIST3(LRI,LRA,LR)
      !WL = WL+2*VINT_CI(LIST+1)-VINT_CI(LIST)       !310:NEOC=2,COE
      WL = 2*VLOP0
      call TRANS_IJKL_INTPOS(LRA,LRI,LR,LR,NXO)
      call PRODAB_2(3,JPEL,IWDL,IWDR,JWL,JWR,WL,JPER,NXO)
      WL = -VLOP0
      call TRANS_IJKL_INTPOS(LRA,LR,LRI,LR,NXO)
      call PRODAB_2(3,JPEL,IWDL,IWDR,JWL,JWR,WL,JPER,NXO)

    end do
    do LRK=norb_dz+1,LRA
      !LIST = LIST3(LRI,LRA,LRK)
      KCOE = LPCOE(LRK)
      call NEOC(KCOE,NOCC,TCOE)
      !WL = WL+NOCC*(VINT_CI(LIST+1)+TCOE*VINT_CI(LIST))
      WL = VLOP0*NOCC
      call TRANS_IJKL_INTPOS(LRA,LRI,LRK,LRK,NXO)
      call PRODAB_2(3,JPEL,IWDL,IWDR,JWL,JWR,WL,JPER,NXO)
      WL = VLOP0*NOCC*TCOE
      call TRANS_IJKL_INTPOS(LRA,LRK,LRI,LRK,NXO)
      call PRODAB_2(3,JPEL,IWDL,IWDR,JWL,JWR,WL,JPER,NXO)

    end do
    !WL = WL*VLOP0
    ! T1D1(15-5) D&rl(33)B^l(13)C'(11)
    VLOP0 = -W0*W0TD5
    do LRK=1,LRI-1
      !LIST = LIST3(LRI,LRA,LRK)
      !WL = WL-VLOP0*(2*VINT_CI(LIST+1)-vint_ci(LIST))
      WL = -2*VLOP0
      call TRANS_IJKL_INTPOS(LRA,LRI,LRK,LRK,NXO)
      call PRODAB_2(3,JPEL,IWDL,IWDR,JWL,JWR,WL,JPER,NXO)
      WL = VLOP0
      call TRANS_IJKL_INTPOS(LRA,LRK,LRI,LRK,NXO)
      call PRODAB_2(3,JPEL,IWDL,IWDR,JWL,JWR,WL,JPER,NXO)

    end do
    !call PRODAB(3,JPEL,IWDL,IWDR,JWL,JWR,WL,JPER)
  end do
  !---------------------------------------------------------------------
  do LRD=NORB_FRZ+1,LRI-1
    LMD = LSM_INN(LRD)
    if (LMD /= JMR) cycle
    IWDL = JUST(LRD,LRI)
    IWDR = JUD(LRD)
    ! T1D1(15-1) (11)A&(13)
    VLOP0 = W0*W0TD1
    !LIST = LIST3(LRI,LRA,LRI)
    !WL = VLOP0*(VOINT(LRI,LRA)+VINT_CI(LIST))             !310,act_
    WL = VLOP0
    call PRODAB_1(3,JPEL,IWDL,IWDR,JWL,JWR,WL,JPER,LRI,LRA)
    call TRANS_IJKL_INTPOS(LRA,LRI,LRI,LRI,NXO)
    call PRODAB_2(3,JPEL,IWDL,IWDR,JWL,JWR,WL,JPER,NXO)

    do LR=LRI+1,NORB_DZ
      !LIST = LIST3(LRI,LRA,LR)
      !WL = WL+VLOP0*(2*VINT_CI(LIST+1)-VINT_CI(LIST)) !  310:NEOC=2
      WL = 2*VLOP0
      call TRANS_IJKL_INTPOS(LRA,LRI,LR,LR,NXO)
      call PRODAB_2(3,JPEL,IWDL,IWDR,JWL,JWR,WL,JPER,NXO)
      WL = -VLOP0
      call TRANS_IJKL_INTPOS(LRA,LR,LRI,LR,NXO)
      call PRODAB_2(3,JPEL,IWDL,IWDR,JWL,JWR,WL,JPER,NXO)

    end do
    do LRK=norb_dz+1,LRA
      !LIST = LIST3(LRI,LRA,LRK)
      KCOE = LPCOE(LRK)
      call NEOC(KCOE,NOCC,TCOE)
      !WL = WL+VLOP0*NOCC*(VINT_CI(LIST+1)+TCOE*VINT_CI(LIST))
      WL = VLOP0*NOCC
      call TRANS_IJKL_INTPOS(LRA,LRI,LRK,LRK,NXO)
      call PRODAB_2(3,JPEL,IWDL,IWDR,JWL,JWR,WL,JPER,NXO)
      WL = VLOP0*NOCC*TCOE
      call TRANS_IJKL_INTPOS(LRA,LRK,LRI,LRK,NXO)
      call PRODAB_2(3,JPEL,IWDL,IWDR,JWL,JWR,WL,JPER,NXO)

    end do
    !WL = WL*VLOP0
    ! T1D1(15-4) D&r&l(11)B^l(13)
    VLOP0 = W0*W0TD4
    VLOP1 = W1*W0TD4
    !LIST = LIST3(LRI,LRA,LRD)
    !WL = WL+(VLOP0-VLOP1)*vint_ci(LIST)-2*VLOP0*VINT_CI(LIST+1)
    WL = VLOP0-VLOP1
    call TRANS_IJKL_INTPOS(LRA,LRD,LRI,LRD,NXO)
    call PRODAB_2(3,JPEL,IWDL,IWDR,JWL,JWR,WL,JPER,NXO)
    WL = -2*VLOP0
    call TRANS_IJKL_INTPOS(LRA,LRI,LRD,LRD,NXO)
    call PRODAB_2(3,JPEL,IWDL,IWDR,JWL,JWR,WL,JPER,NXO)

    ! T1D1(15-5) D&rl(33)C"(11)B^l(13)
    VLOP0 = W0*W0TD5
    do LRK=1,LRI-1
      if (LRK == LRD) cycle
      !LIST = LIST3(LRI,LRA,LRK)
      !WL = WL+VLOP0*(vint_ci(LIST)-2*VINT_CI(LIST+1))      !4.3
      WL = -2*VLOP0
      call TRANS_IJKL_INTPOS(LRA,LRI,LRK,LRK,NXO)
      call PRODAB_2(3,JPEL,IWDL,IWDR,JWL,JWR,WL,JPER,NXO)
      WL = VLOP0
      call TRANS_IJKL_INTPOS(LRA,LRK,LRI,LRK,NXO)
      call PRODAB_2(3,JPEL,IWDL,IWDR,JWL,JWR,WL,JPER,NXO)

    end do
    !call PRODAB(3,JPEL,IWDL,IWDR,JWL,JWR,WL,JPER)
  end do
end do
do LRI=NORB_FRZ+1,NORB_DZ-1
  LMI = LSM_INN(LRI)
  do LRJ=LRI+1,NORB_DZ
    LMJ = LSM_INN(LRJ)
    LMIJ = Mul(LMI,LMJ)
    if (LMIJ /= JML) cycle
    IWDL = JUST(LRI,LRJ)

    ! T1D1(15-2) A&(13)B&r(13)B^r(31)
    do LRD=LRJ+1,NORB_DZ
      LMD = LSM_INN(LRD)
      if (LMD /= JMR) cycle
      W0TD2 = W0_T1D1(2)
      W1TD2 = W1_T1D1(2)
      NI = mod(LRJ-LRI+NORB_DZ-LRD,2)
      if (NI == 0) W0TD2 = -W0TD2
      if (NI == 0) W1TD2 = -W1TD2

      IWDR = JUD(LRD)
      VLOP0 = W0*W0TD2
      VLOP1 = W1*W1TD2
      !LIST = LIST4(LRI,LRJ,LRD,LRA)
      !WL = VLOP0*(vint_ci(LIST+2)+vint_ci(LIST))-VLOP1*(vint_ci(LIST+2)-vint_ci(LIST))  !1.3
      WL = VLOP0-VLOP1
      call TRANS_IJKL_INTPOS(LRA,LRI,LRJ,LRD,NXO)
      call PRODAB_2(3,JPEL,IWDL,IWDR,JWL,JWR,WL,JPER,NXO)
      WL = VLOP0+VLOP1
      call TRANS_IJKL_INTPOS(LRA,LRJ,LRD,LRI,NXO)
      call PRODAB_2(3,JPEL,IWDL,IWDR,JWL,JWR,WL,JPER,NXO)

      !call PRODAB(3,JPEL,IWDL,IWDR,JWL,JWR,WL,JPER)
    end do
    ! T1D1(15-3) A&(13)B&l(31)B^l(13)
    do LRD=LRI+1,LRJ-1
      LMD = LSM_INN(LRD)
      if (LMD /= JMR) cycle
      IWDR = JUD(LRD)
      W0TD3 = W0_T1D1(3)
      W1TD3 = W1_T1D1(3)
      NI = mod(LRD-LRI+NORB_DZ-LRJ,2)
      if (NI == 0) W0TD3 = -W0TD3
      if (NI == 0) W1TD3 = -W1TD3
      VLOP0 = W0*W0TD3                !D6-8
      VLOP1 = W1*W1TD3
      !LIST = LIST4(LRI,LRD,LRJ,LRA)
      !WL = VLOP0*(vint_ci(LIST+2)-2*vint_ci(LIST+1))-VLOP1*vint_ci(LIST+2)      !1.2
      WL = VLOP0-VLOP1
      call TRANS_IJKL_INTPOS(LRA,LRI,LRD,LRJ,NXO)
      call PRODAB_2(3,JPEL,IWDL,IWDR,JWL,JWR,WL,JPER,NXO)
      WL = -2*VLOP0
      call TRANS_IJKL_INTPOS(LRA,LRJ,LRD,LRI,NXO)
      call PRODAB_2(3,JPEL,IWDL,IWDR,JWL,JWR,WL,JPER,NXO)

      !call PRODAB(3,JPEL,IWDL,IWDR,JWL,JWR,WL,JPER)
    end do
  end do
end do

return

end subroutine TTDD_HEAD_DBL_TAIL_ACT_G
