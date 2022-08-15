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

subroutine ss_head_dbl_tail_act(lra)

use gugaci_global, only: jb_sys, jml, jmr, jpad, jpadl, jpel, jper, just, jwl, jwr, line, lrs, lsm_inn, norb_dz, norb_frz, w0, &
                         w0_ss, w1, w1_ss
use Symmetry_Info, only: Mul
use Constants, only: Zero
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: lra
integer(kind=iwp) :: iwdl, iwdr, jmlr, lmi, lmij, lmj, lmk, lmki, lmkj, lr0, lri, lrj, lrk, ni
real(kind=wp) :: vlop0, vlop1, w0ss1, w0ss10, w0ss11, w0ss12, w0ss13, w0ss14, w0ss15, w0ss16, w0ss17, w0ss18, w0ss2, w0ss20, &
                 w0ss3, w0ss4, w0ss5, w0ss6, w0ss7, w0ss8, w0ss9, w1ss1, w1ss10, w1ss11, w1ss12, w1ss13, w1ss14, w1ss15, w1ss16, &
                 w1ss17, w1ss18, w1ss2, w1ss3, w1ss4, w1ss5, w1ss6, w1ss7, w1ss8, w1ss9, wl, wl1, wl2, wl3, wl4

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

! ss(1-1)  ar(01)-bl(32)-
! ss(1-2)  ar(02)-bl(31)-
! ss(1-3)  ar(13)-bl(20)-
! ss(1-4)  ar(23)-bl(10)-
! ss(1-5)  (22)-ar(13)-bl(31)-
! ss(1-6)  (11)-ar(23)-bl(32)-
! ss(1-7)  ar(13)-c'(21)-bl(32)-
! ss(1-8)  ar(13)-c'(22)-bl(31)-
! ss(1-9)  ar(23)-c'(11)-bl(32)-
! ss(1-10) ar(23)-c'(12)-bl(31)-
! ss(1-11) ar(13)-bl(31)-c"(22)-
! ss(1-12) ar(13)-bl(32)-c"(21)-
! ss(1-13) ar(23)-bl(31)-c"(12)-
! ss(1-14) ar(23)-bl(32)-c"(11)-
! ss(1-15) (22)-drl(11)-
! ss(1-16) (11)-drl(22)-
! ss(1-17) drl(22)-c"(11)-
! ss(1-18) drl(11)-c"(22)-
! ss(1-19) drl(12)-c"(21)-
! ss(1-20) drl(33)-c"(00)-
! ss(1-20) drl(33)-c"(11)-c"(22)-
! ss(1-20) (11)drl(33)-c"(22)-
! ss(1-20) (11)(22)drl(33)-
! ss(1-20) drl(33)-c"(22)-c"(11)-
! ss(1-20) (22)drl(33)-c"(11)-
! ss(1-20) (22)(11)drl(33)-
! ss(1-20) (21)(12)drl(33)-
do lri=norb_frz+1,norb_dz-1
  lmi = lsm_inn(lri)
  do lrj=lri+1,norb_dz
    lmj = lsm_inn(lrj)
    lmij = Mul(lmi,lmj)
    jmlr = Mul(jml,jmr)
    if (lmij /= jmlr) cycle
    w0ss2 = w0_ss(2)
    w1ss2 = w1_ss(2)
    w0ss4 = w0_ss(4)
    w1ss4 = w1_ss(4)
    w0ss5 = w0_ss(5)
    w1ss5 = w1_ss(5)
    w0ss10 = -w0_ss(10)
    w1ss10 = -w1_ss(10)
    w0ss14 = w0_ss(14)
    w1ss14 = w1_ss(14)
    if (jb_sys > 0) then
      w0ss1 = w0_ss(1)
      w1ss1 = w1_ss(1)
      w0ss3 = w0_ss(3)
      w1ss3 = w1_ss(3)
      w0ss6 = w0_ss(6)
      w1ss6 = w1_ss(6)
      w0ss7 = w0_ss(7)
      w1ss7 = w1_ss(7)
      w0ss8 = w0_ss(8)
      w1ss8 = w1_ss(8)
      w0ss9 = w0_ss(9)
      w1ss9 = w1_ss(9)
      w0ss11 = w0_ss(11)
      w1ss11 = w1_ss(11)
      w0ss12 = w0_ss(12)
      w1ss12 = w1_ss(12)
      w0ss13 = w0_ss(13)
      w1ss13 = w1_ss(13)
    end if
    ni = mod(lrj-lri,2)
    if (ni == 0) then
      w0ss2 = -w0ss2
      w1ss2 = -w1ss2
      w0ss4 = -w0ss4
      w1ss4 = -w1ss4
      w0ss5 = -w0ss5
      w1ss5 = -w1ss5
      w0ss10 = -w0ss10
      w1ss10 = -w1ss10
      w0ss14 = -w0ss14
      w1ss14 = -w1ss14
      if (jb_sys > 0) then
        w0ss1 = -w0ss1
        w1ss1 = -w1ss1
        w0ss3 = -w0ss3
        w1ss3 = -w1ss3
        w0ss6 = -w0ss6
        w1ss6 = -w1ss6
        w0ss7 = -w0ss7
        w1ss7 = -w1ss7
        w0ss8 = -w0ss8
        w1ss8 = -w1ss8
        w0ss9 = -w0ss9
        w1ss9 = -w1ss9
        w0ss11 = -w0ss11
        w1ss11 = -w1ss11
        w0ss12 = -w0ss12
        w1ss12 = -w1ss12
        w0ss13 = -w0ss13
        w1ss13 = -w1ss13
      end if
    end if

    if ((jml == 1) .and. (lmij == jmr)) then
      iwdl = just(lri,lri)
      ! ss(1-1)   ar(01)-bl(32)-
      if (jb_sys > 0) then
        iwdr = just(lrj,lri)
        vlop0 = w0*w0ss1
        vlop1 = w1*w1ss1
        if (line == 26) then    !lri,lrj,lra
          call comp_loop(5,lri,lrj,lrs,lra,vlop0,vlop1,wl)
        end if
        if (line == 28) then    !lri,lrj,lrs,lra
          call comp_loop(7,lri,lrj,lrs,lra,vlop0,vlop1,wl)
        end if
        if (line == 29) then    !lri,lrj,lrs,lra
          call comp_loop(6,lri,lrj,lrs,lra,vlop0,vlop1,wl)
        end if
        call prodab(3,jpel,iwdl,iwdr,jwl,jwr,wl,jper)
      end if
      ! ss(1-2)  ar(02)-bl(31)-
      iwdr = just(lri,lrj)
      vlop0 = w0*w0ss2
      vlop1 = w1*w1ss2
      if (line == 26) then    !lri,lrj,lra
        call comp_loop(5,lri,lrj,lrs,lra,vlop0,vlop1,wl)
      end if
      if (line == 28) then    !lri,lrj,lrs,lra
        call comp_loop(7,lri,lrj,lrs,lra,vlop0,vlop1,wl)
      end if
      if (line == 29) then    !lri,lrj,lrs,lra
        call comp_loop(6,lri,lrj,lrs,lra,vlop0,vlop1,wl)
      end if
      call prodab(3,jpel,iwdl,iwdr,jwl,jwr,wl,jper)
    end if

    if ((jmr == 1) .and. (lmij == jml)) then
      iwdr = just(lrj,lrj)
      ! ss(1-3)  ar(13)-bl(20)-
      if (jb_sys > 0) then
        iwdl = just(lrj,lri)
        vlop0 = w0*w0ss3
        vlop1 = w1*w1ss3
        if (line == 26) then    !lri,lrj,lra
          call comp_loop(5,lri,lrj,lrs,lra,vlop0,vlop1,wl)
        end if
        if (line == 28) then    !lri,lrj,lrs,lra
          call comp_loop(7,lri,lrj,lrs,lra,vlop0,vlop1,wl)
        end if
        if (line == 29) then    !lri,lrj,lrs,lra
          call comp_loop(6,lri,lrj,lrs,lra,vlop0,vlop1,wl)
        end if
        call prodab(3,jpel,iwdl,iwdr,jwl,jwr,wl,jper)
      end if
      ! ss(1-4)  ar(23)-bl(10)-        act -c"-                  ! iprad
      iwdl = just(lri,lrj)
      vlop0 = w0*w0ss4
      vlop1 = w1*w1ss4
      if (line == 26) then    !lri,lrj,lra
        call comp_loop(5,lri,lrj,lrs,lra,vlop0,vlop1,wl)
      end if
      if (line == 28) then    !lri,lrj,lrs,lra
        call comp_loop(7,lri,lrj,lrs,lra,vlop0,vlop1,wl)
      end if
      if (line == 29) then    !lri,lrj,lrs,lra
        call comp_loop(6,lri,lrj,lrs,lra,vlop0,vlop1,wl)
      end if
      call prodab(3,jpel,iwdl,iwdr,jwl,jwr,wl,jper)
    end if

    vlop0 = w0*w0ss5
    vlop1 = w1*w1ss5
    if (line == 26) then    !lri,lrj,lra
      call comp_loop(5,lri,lrj,lrs,lra,vlop0,vlop1,wl)
    end if
    if (line == 28) then    !lri,lrj,lrs,lra
      call comp_loop(7,lri,lrj,lrs,lra,vlop0,vlop1,wl)
    end if
    if (line == 29) then    !lri,lrj,lrs,lra
      call comp_loop(6,lri,lrj,lrs,lra,vlop0,vlop1,wl)
    end if
    if (jb_sys > 0) then
      vlop0 = w0*w0ss6
      vlop1 = w1*w1ss6
      if (line == 26) then    !lri,lrj,lra
        call comp_loop(5,lri,lrj,lrs,lra,vlop0,vlop1,wl1)
      end if
      if (line == 28) then    !lri,lrj,lrs,lra
        call comp_loop(7,lri,lrj,lrs,lra,vlop0,vlop1,wl1)
      end if
      if (line == 29) then    !lri,lrj,lrs,lra
        call comp_loop(6,lri,lrj,lrs,lra,vlop0,vlop1,wl1)
      end if
    end if
    do lrk=norb_frz+1,lri-1
      lmk = lsm_inn(lrk)
      lmki = Mul(lmk,lmi)
      lmkj = Mul(lmk,lmj)
      if ((lmki == jml) .and. (lmkj == jmr)) then
        ! ss(1-5)  (22)-ar(13)-bl(31)-
        iwdl = just(lrk,lri)
        iwdr = just(lrk,lrj)
        call prodab(3,jpel,iwdl,iwdr,jwl,jwr,wl,jper)
        if (jb_sys > 0) then
          ! ss(1-6)  (11)-ar(23)-bl(32)-
          iwdl = just(lri,lrk)
          iwdr = just(lrj,lrk)
          call prodab(3,jpel,iwdl,iwdr,jwl,jwr,wl1,jper)
        end if
      end if
    end do

    if (jb_sys > 0) then
      vlop0 = w0*w0ss7
      vlop1 = w1*w1ss7
      if (line == 26) then    !lri,lrj,lra
        call comp_loop(5,lri,lrj,lrs,lra,vlop0,vlop1,wl1)
      end if
      if (line == 28) then    !lri,lrj,lrs,lra
        call comp_loop(7,lri,lrj,lrs,lra,vlop0,vlop1,wl1)
      end if
      if (line == 29) then    !lri,lrj,lrs,lra
        call comp_loop(6,lri,lrj,lrs,lra,vlop0,vlop1,wl1)
      end if
      vlop0 = w0*w0ss8
      vlop1 = w1*w1ss8
      if (line == 26) then    !lri,lrj,lra
        call comp_loop(5,lri,lrj,lrs,lra,vlop0,vlop1,wl2)
      end if
      if (line == 28) then    !lri,lrj,lrs,lra
        call comp_loop(7,lri,lrj,lrs,lra,vlop0,vlop1,wl2)
      end if
      if (line == 29) then    !lri,lrj,lrs,lra
        call comp_loop(6,lri,lrj,lrs,lra,vlop0,vlop1,wl2)
      end if
      vlop0 = w0*w0ss9
      vlop1 = w1*w1ss9
      if (line == 26) then    !lri,lrj,lra
        call comp_loop(5,lri,lrj,lrs,lra,vlop0,vlop1,wl3)
      end if
      if (line == 28) then    !lri,lrj,lrs,lra
        call comp_loop(7,lri,lrj,lrs,lra,vlop0,vlop1,wl3)
      end if
      if (line == 29) then    !lri,lrj,lrs,lra
        call comp_loop(6,lri,lrj,lrs,lra,vlop0,vlop1,wl3)
      end if
    end if
    vlop0 = w0*w0ss10
    vlop1 = w1*w1ss10
    if (line == 26) then    !lri,lrj,lra
      call comp_loop(5,lri,lrj,lrs,lra,vlop0,vlop1,wl4)
    end if
    if (line == 28) then    !lri,lrj,lrs,lra
      call comp_loop(7,lri,lrj,lrs,lra,vlop0,vlop1,wl4)
    end if
    if (line == 29) then    !lri,lrj,lrs,lra
      call comp_loop(6,lri,lrj,lrs,lra,vlop0,vlop1,wl4)
    end if

    do lrk=lri+1,lrj-1
      lmk = lsm_inn(lrk)
      lmki = Mul(lmk,lmi)
      lmkj = Mul(lmk,lmj)
      if ((lmki == jml) .and. (lmkj == jmr)) then
        if (jb_sys > 0) then
          ! ss(1-7)  ar(13)-c'(21)-bl(32)-
          iwdl = just(lrk,lri)
          iwdr = just(lrj,lrk)
          call prodab(3,jpel,iwdl,iwdr,jwl,jwr,-wl1,jper)
          ! ss(1-8)  ar(13)-c'(22)-bl(31)-
          iwdl = just(lrk,lri)
          iwdr = just(lrk,lrj)
          call prodab(3,jpel,iwdl,iwdr,jwl,jwr,-wl2,jper)
          ! ss(1-9)  ar(23)-c'(11)-bl(32)-
          iwdl = just(lri,lrk)
          iwdr = just(lrj,lrk)
          !wl = (vlop0-vlop1)*vint_ci(list)-2*vlop0*vint_ci(list+1)
          call prodab(3,jpel,iwdl,iwdr,jwl,jwr,-wl3,jper)
        end if
        ! ss(1-10) ar(23)-c'(12)-bl(31)-
        iwdl = just(lri,lrk)
        iwdr = just(lrk,lrj)
        call prodab(3,jpel,iwdl,iwdr,jwl,jwr,wl4,jper)
      end if
    end do

    if (jb_sys > 0) then
      vlop0 = w0*w0ss11
      vlop1 = w1*w1ss11
      if (line == 26) then    !lri,lrj,lra
        call comp_loop(5,lri,lrj,lrs,lra,vlop0,vlop1,wl1)
      end if
      if (line == 28) then    !lri,lrj,lrs,lra
        call comp_loop(7,lri,lrj,lrs,lra,vlop0,vlop1,wl1)
      end if
      if (line == 29) then    !lri,lrj,lrs,lra
        call comp_loop(6,lri,lrj,lrs,lra,vlop0,vlop1,wl1)
      end if
      vlop0 = w0*w0ss12
      vlop1 = w1*w1ss12
      if (line == 26) then    !lri,lrj,lra
        call comp_loop(5,lri,lrj,lrs,lra,vlop0,vlop1,wl2)
      end if
      if (line == 28) then    !lri,lrj,lrs,lra
        call comp_loop(7,lri,lrj,lrs,lra,vlop0,vlop1,wl2)
      end if
      if (line == 29) then    !lri,lrj,lrs,lra
        call comp_loop(6,lri,lrj,lrs,lra,vlop0,vlop1,wl2)
      end if
      vlop0 = w0*w0ss13
      vlop1 = w1*w1ss13
      if (line == 26) then    !lri,lrj,lra
        call comp_loop(5,lri,lrj,lrs,lra,vlop0,vlop1,wl3)
      end if
      if (line == 28) then    !lri,lrj,lrs,lra
        call comp_loop(7,lri,lrj,lrs,lra,vlop0,vlop1,wl3)
      end if
      if (line == 29) then    !lri,lrj,lrs,lra
        call comp_loop(6,lri,lrj,lrs,lra,vlop0,vlop1,wl3)
      end if
    end if
    vlop0 = w0*w0ss14
    vlop1 = w1*w1ss14
    if (line == 26) then    !lri,lrj,lra
      call comp_loop(5,lri,lrj,lrs,lra,vlop0,vlop1,wl4)
    end if
    if (line == 28) then    !lri,lrj,lrs,lra
      call comp_loop(7,lri,lrj,lrs,lra,vlop0,vlop1,wl4)
    end if
    if (line == 29) then    !lri,lrj,lrs,lra
      call comp_loop(6,lri,lrj,lrs,lra,vlop0,vlop1,wl4)
    end if

    do lrk=lrj+1,norb_dz
      lmk = lsm_inn(lrk)
      lmki = Mul(lmk,lmi)
      lmkj = Mul(lmk,lmj)
      if ((lmki == jml) .and. (lmkj == jmr)) then
        if (jb_sys > 0) then
          ! ss(1-11) ar(13)-bl(31)-c"(22)-
          iwdl = just(lrk,lri)
          iwdr = just(lrk,lrj)
          call prodab(3,jpel,iwdl,iwdr,jwl,jwr,wl1,jper)
          ! ss(1-12) ar(13)-bl(32)-c"(21)-
          iwdl = just(lrk,lri)
          iwdr = just(lrj,lrk)
          call prodab(3,jpel,iwdl,iwdr,jwl,jwr,wl2,jper)
          ! ss(1-13) ar(23)-bl(31)-c"(12)-
          iwdl = just(lri,lrk)
          iwdr = just(lrk,lrj)
          call prodab(3,jpel,iwdl,iwdr,jwl,jwr,wl3,jper)
        end if
        ! ss(1-14) ar(23)-bl(32)-c"(11)- act -c"-
        iwdl = just(lri,lrk)
        iwdr = just(lrj,lrk)
        call prodab(3,jpel,iwdl,iwdr,jwl,jwr,wl4,jper)
      end if
    end do
  end do
end do

if (jpad /= jpadl) return

!if ((jb_sys > 0) .or. (jwl >= jwr)) then
if (jb_sys > 0) then
  if ((jpad > 17) .and. (jpad < 25)) then
    do lri=norb_frz+1,norb_dz
      lmi = lsm_inn(lri)
      do lrj=lri+1,norb_dz
        lmj = lsm_inn(lrj)
        lmij = Mul(lmi,lmj)
        if ((lmij /= jml) .or. (jml /= jmr)) cycle
        iwdl = just(lrj,lri)
        iwdr = just(lri,lrj)
        ! ss(1-19) drl(12)-c"(21)-
        vlop0 = w0*w0_ss(19)
        vlop1 = w1*w1_ss(19)
        if (line == 26) then    !lri,lra
          call comp_loop(9,lri,0,0,lra,vlop0,vlop1,wl)
        end if
        if (line == 28) then    !lri,lrs,lra
          call comp_loop(12,lri,0,lrs,lra,vlop0,vlop1,wl)
        end if
        if (line == 29) then    !lri,lrs,lra
          call comp_loop(11,lri,0,lrs,lra,vlop0,vlop1,wl)
        end if
        !wl = (vlop0-vlop1)*voint(lri,lrb)
        call prodab(3,jpel,iwdl,iwdr,jwl,jwr,wl,jper)
      end do
    end do
  end if
end if
if (jwl >= jwr) return

w0ss15 = w0_ss(15)
w1ss15 = w1_ss(15)
w0ss17 = w0_ss(17)
w1ss17 = w1_ss(17)
w0ss20 = w0_ss(20)
if (jb_sys > 0) then
  w0ss16 = w0_ss(16)
  w1ss16 = w1_ss(16)
  w0ss18 = w0_ss(18)
  w1ss18 = w1_ss(18)
end if

if ((jml == 1) .and. (jmr == 1)) then
  ! ss(1-20) drl(33)-c"(00)-                          ! ipl(r)ad=1 or =ns
  do lr0=norb_frz+1,norb_dz
    iwdl = just(lr0,lr0)
    iwdr = iwdl
    vlop0 = w0*w0_ss(20)
    vlop1 = Zero
    wl = Zero
    do lrk=1,norb_dz
      if (lrk == lr0) cycle
      if (line == 26) then    !lrk,lra
        call comp_loop(9,lrk,0,0,lra,vlop0,vlop1,wl1)
      end if
      if (line == 28) then    !lrk,lrs,lra    lrs,lra,lrk
        !call comp_loop(12,lrk,lrs,lra,lra,vlop0,vlop1,wl1)
        call comp_loop(12,lrk,0,lrs,lra,vlop0,vlop1,wl1)
      end if
      if (line == 29) then    !lrk,lrs,lra
        call comp_loop(11,lrk,0,lrs,lra,vlop0,vlop1,wl1)
        !call comp_loop(11,lrk,lrs,lra,lra,vlop0,vlop1,wl1)
      end if
      wl = wl+wl1
    end do
    call prodab(3,jpel,iwdl,iwdr,jwl,jwr,wl,jper)
  end do
end if
do lri=norb_frz+1,norb_dz
  lmi = lsm_inn(lri)
  do lrj=lri+1,norb_dz
    lmj = lsm_inn(lrj)
    lmij = Mul(lmi,lmj)
    if ((lmij /= jml) .or. (jml /= jmr)) cycle
    if ((jwl >= jwr) .and. (jb_sys == 0)) cycle
    wl = Zero
    iwdl = just(lri,lrj)
    iwdr = iwdl
    ! ss(1-15) (22)-drl(11)-
    vlop0 = w0*w0ss15
    vlop1 = w1*w1ss15
    wl = Zero
    if (line == 26) then    !lrj,lra
      call comp_loop(9,lrj,0,0,lra,vlop0,vlop1,wl1)
    end if
    if (line == 28) then    !lrj,lrs,lra
      call comp_loop(12,lrj,0,lrs,lra,vlop0,vlop1,wl1)
    end if
    if (line == 29) then    !lrj,lrs,lra
      call comp_loop(11,lrj,0,lrs,lra,vlop0,vlop1,wl1)
    end if
    wl = wl+wl1
    !wl = wl+(vlop0-vlop1)*voint(lrj,lrb)
    ! ss(1-17) drl(22)-c"(11)-
    vlop0 = w0*w0ss17
    vlop1 = w1*w1ss17
    if (line == 26) then    !lri,lra
      call comp_loop(9,lri,0,0,lra,vlop0,vlop1,wl1)
    end if
    if (line == 28) then    !lri,lrs,lra
      call comp_loop(12,lri,0,lrs,lra,vlop0,vlop1,wl1)
    end if
    if (line == 29) then    !lri,lrs,lra
      call comp_loop(11,lri,0,lrs,lra,vlop0,vlop1,wl1)
    end if
    wl = wl+wl1
    !wl = wl+(vlop0-vlop1)*voint(lri,lrb)
    ! ss(1-20) (22)(11)drl(33)-
    ! ss(1-20) (22)drl(33)-c"(11)-
    ! ss(1-20) drl(33)-c"(22)-c"(11)-
    vlop0 = w0*w0ss20
    vlop1 = Zero
    do lrk=1,norb_dz
      if (lrk == lri) cycle
      if (lrk == lrj) cycle
      if (line == 26) then    !lrk,lra
        call comp_loop(9,lrk,0,0,lra,vlop0,vlop1,wl1)
      end if
      if (line == 28) then    !lrk,lrs,lra
        call comp_loop(12,lrk,0,lrs,lra,vlop0,vlop1,wl1)
      end if
      if (line == 29) then    !lrk,lrs,lra
        call comp_loop(11,lrk,0,lrs,lra,vlop0,vlop1,wl1)
      end if
      wl = wl+wl1
      !wl = wl+vlop0*voint(lrk,lrb)
    end do
    call prodab(3,jpel,iwdl,iwdr,jwl,jwr,wl,jper)
    if (jb_sys > 0) then
      iwdl = just(lrj,lri)
      iwdr = iwdl
      !ss(1-16) (11)-drl(22)-
      vlop0 = w0*w0ss16
      vlop1 = w1*w1ss16
      wl = Zero
      if (line == 26) then    !lrj,lra
        call comp_loop(9,lrj,0,0,lra,vlop0,vlop1,wl1)
      end if
      if (line == 28) then    !lrj,lrs,lra
        call comp_loop(12,lrj,0,lrs,lra,vlop0,vlop1,wl1)
      end if
      if (line == 29) then    !lrj,lrs,lra
        call comp_loop(11,lrj,0,lrs,lra,vlop0,vlop1,wl1)
      end if
      wl = wl+wl1
      !wl = (vlop0-vlop1)*voint(lrj,lrb)
      ! ss(1-18) drl(11)-c"(22)-
      vlop0 = w0*w0ss18
      vlop1 = w1*w1ss18
      if (line == 26) then    !lri,lra
        call comp_loop(9,lri,0,0,lra,vlop0,vlop1,wl1)
      end if
      if (line == 28) then    !lri,lrs,lra
        call comp_loop(12,lri,0,lrs,lra,vlop0,vlop1,wl1)
      end if
      if (line == 29) then    !lri,lrs,lra
        call comp_loop(11,lri,0,lrs,lra,vlop0,vlop1,wl1)
      end if
      wl = wl+wl1
      !wl = wl+(vlop0-vlop1)*voint(lri,lrb)
      ! ss(1-20) drl(33)-c"(11)-c"(22)-
      ! ss(1-20) (11)drl(33)-c"(22)-
      ! ss(1-20) (11)(22)drl(33)-
      vlop0 = w0*w0ss20
      vlop1 = Zero
      do lrk=1,norb_dz
        if (lrk == lri) cycle
        if (lrk == lrj) cycle
        if (line == 26) then    !lrk,lra
          call comp_loop(9,lrk,0,0,lra,vlop0,vlop1,wl1)
        end if
        if (line == 28) then    !lrk,lrs,lra
          call comp_loop(12,lrk,0,lrs,lra,vlop0,vlop1,wl1)
        end if
        if (line == 29) then    !lrk,lrs,lra
          call comp_loop(11,lrk,0,lrs,lra,vlop0,vlop1,wl1)
        end if
        wl = wl+wl1
        !wl = wl+vlop0*voint(lrk,lrb)
      end do
      call prodab(3,jpel,iwdl,iwdr,jwl,jwr,wl,jper)
    end if
  end do
end do

return

end subroutine ss_head_dbl_tail_act

subroutine st_head_dbl_tail_act(lra)

use gugaci_global, only: jb_sys, jml, jmr, jpel, jper, just, jwl, jwr, line, lrs, lsm_inn, norb_dz, norb_frz, w1, w1_st
use Symmetry_Info, only: Mul
use Constants, only: Zero
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: lra
integer(kind=iwp) :: iwds, iwdt, jmlr, lmi, lmij, lmj, lmk, lri, lrj, lrk, ni
real(kind=wp) :: vlop0, vlop1, w1st1, w1st2, w1st3, w1st4, wl, wl1

! st(2-1) ar(02)-bl(32)-
! st(2-2) (22)ar(13)-bl(32)-
! st(2-4) ar(23)-c'(12)-bl(32)-
! st(2-4) ar(23)-bl(32)-c'(12)-
! st(2-5) (22)drl(12)-
! st(2-6) drl(22)-c"(12)-
do lri=norb_frz+1,norb_dz
  lmi = lsm_inn(lri)
  do lrj=lri+1,norb_dz
    lmj = lsm_inn(lrj)
    lmij = Mul(lmi,lmj)
    lmij = Mul(lmij,1)
    if ((jml == jmr) .and. (lmij == jml)) then
      iwds = just(lri,lrj)
      iwdt = iwds
      ! st(2-5) (22)drl(12)-
      vlop1 = w1*w1_st(5)             !d2-5
      vlop0 = Zero
      if (line == 26) then    !lrj,lra
        call comp_loop(9,lrj,0,0,lra,vlop0,vlop1,wl)
      end if
      if (line == 28) then    !lrj,lrs,lra
        call comp_loop(12,lrj,0,lrs,lra,vlop0,vlop1,wl)
      end if
      if (line == 29) then    !lrj,lrs,lra
        call comp_loop(11,lrj,0,lrs,lra,vlop0,vlop1,wl)
      end if
      ! st(2-6) drl(22)-c"(12)-
      vlop1 = w1*w1_st(6)             !d2-6
      vlop0 = Zero
      if (line == 26) then    !lri,lra
        call comp_loop(9,lri,0,0,lra,vlop0,vlop1,wl1)
      end if
      if (line == 28) then    !lri,lrs,lra
        call comp_loop(12,lri,0,lrs,lra,vlop0,vlop1,wl1)
      end if
      if (line == 29) then    !lri,lrs,lra
        call comp_loop(11,lri,0,lrs,lra,vlop0,vlop1,wl1)
      end if
      wl = wl+wl1
      !list = list3(lra,lrb,lri)
      !wl = wl-vlop1*vint_ci(list)    !4.3 vlop0=0        !!!!!
      call prodab(3,jpel,iwds,iwdt,jwl,jwr,wl,jper)
      ! st(2-7) drl(12)-c"(22)-
      if (jb_sys > 0) then
        iwds = just(lrj,lri)
        iwdt = just(lri,lrj)
        vlop1 = w1*w1_st(7)
        vlop0 = Zero             !d2-6
        !list = list3(lra,lrb,lri)
        !wl = wl-vlop1*vint_ci(list)    !4.3 vlop0=0        !!!!!
        if (line == 26) then    !lri,lra
          call comp_loop(9,lri,0,0,lra,vlop0,vlop1,wl)
        end if
        if (line == 28) then    !lri,lrs,lra
          call comp_loop(12,lri,0,lrs,lra,vlop0,vlop1,wl)
        end if
        if (line == 29) then    !lri,lrs,lra
          call comp_loop(11,lri,0,lrs,lra,vlop0,vlop1,wl)
        end if
        call prodab(3,jpel,iwds,iwdt,jwl,jwr,wl,jper)
      end if
    end if

    jmlr = Mul(jml,jmr)
    if (lmij /= jmlr) cycle
    w1st1 = w1_st(1)
    w1st2 = w1_st(2)
    w1st3 = w1_st(3)
    w1st4 = w1_st(4)
    ni = mod(lrj-lri,2)
    if (ni == 0) then
      w1st1 = -w1st1
      w1st2 = -w1st2
      w1st3 = -w1st3
      w1st4 = -w1st4
    end if
    if (jml == 1) then
      ! st(2-1) ar(02)-bl(32)-
      iwds = just(lri,lri)
      iwdt = just(lri,lrj)
      vlop1 = w1*w1st1
      vlop0 = Zero
      !list = list4(lri,lrj,lra,lrb)
      if (line == 26) then    !lri,lrj,lra
        call comp_loop(5,lri,lrj,lrs,lra,vlop0,vlop1,wl)
      end if
      if (line == 28) then     !lri,lrj,lrs,lra
        call comp_loop(7,lri,lrj,lrs,lra,vlop0,vlop1,wl)
      end if
      if (line == 29) then     !lri,lrj,lrs,lra
        call comp_loop(6,lri,lrj,lrs,lra,vlop0,vlop1,wl)
      end if
      !wl = -vlop1*vint_ci(list)    !1.1 vlop0=0        !!!!!
      call prodab(3,jpel,iwds,iwdt,jwl,jwr,wl,jper)
    end if

    ! st(2-2) (22)ar(13)-bl(32)-
    vlop1 = w1*w1st2
    !wl = -vlop1*vint_ci(list)    !1.1
    vlop0 = Zero
    if (line == 26) then    !lri,lrj,lra
      call comp_loop(5,lri,lrj,lrs,lra,vlop0,vlop1,wl)
    end if
    if (line == 28) then    !lri,lrj,lrs,lra
      call comp_loop(7,lri,lrj,lrs,lra,vlop0,vlop1,wl)
    end if
    if (line == 29) then    !lri,lrj,lrs,lra
      call comp_loop(6,lri,lrj,lrs,lra,vlop0,vlop1,wl)
    end if
    do lrk=norb_frz+1,lri-1
      lmk = lsm_inn(lrk)
      if (Mul(lmk,lmi) /= jml) cycle
      if (Mul(lmk,lmj) /= jmr) cycle
      iwds = just(lrk,lri)
      iwdt = just(lrk,lrj)
      call prodab(3,jpel,iwds,iwdt,jwl,jwr,wl,jper)
    end do
    ! st(2-3) ar(13)-bl(32)-c'(22)-
    ! st(2-4) ar(23)-bl(32)-c'(12)-
    vlop1 = w1*w1st4
    vlop0 = Zero
    if (line == 26) then    !lri,lrj,lra
      call comp_loop(5,lri,lrj,lrs,lra,vlop0,vlop1,wl)
    end if
    if (line == 28) then    !lri,lrj,lrs,lra
      call comp_loop(7,lri,lrj,lrs,lra,vlop0,vlop1,wl)
    end if
    if (line == 29) then    !lri,lrj,lrs,lra
      call comp_loop(6,lri,lrj,lrs,lra,vlop0,vlop1,wl)
    end if
    !wl = -vlop1*vint_ci(list)    !1.1 vlop0=0
    if (jb_sys > 0) then
      vlop1 = w1*w1st3
      vlop0 = Zero
      if (line == 26) then    !lri,lrj,lra
        call comp_loop(5,lri,lrj,lrs,lra,vlop0,vlop1,wl1)
      end if
      if (line == 28) then    !lri,lrj,lrs,lra
        call comp_loop(7,lri,lrj,lrs,lra,vlop0,vlop1,wl1)
      end if
      if (line == 29) then    !lri,lrj,lrs,lra
        call comp_loop(6,lri,lrj,lrs,lra,vlop0,vlop1,wl1)
      end if
    end if
    !wl1 = -vlop1*vint_ci(list)
    do lrk=lrj+1,norb_dz
      lmk = lsm_inn(lrk)
      if (Mul(lmk,lmi) /= jml) cycle
      if (Mul(lmk,lmj) /= jmr) cycle
      iwds = just(lri,lrk)
      iwdt = just(lrj,lrk)
      call prodab(3,jpel,iwds,iwdt,jwl,jwr,wl,jper)
      if (jb_sys > 0) then
        iwds = just(lrk,lri)
        call prodab(3,jpel,iwds,iwdt,jwl,jwr,wl1,jper)
      end if
    end do
    ! st(2-4) ar(23)-c'(12)-bl(32)-
    do lrk=lri+1,lrj-1
      lmk = lsm_inn(lrk)
      if (Mul(lmk,lmi) /= jml) cycle
      if (Mul(lmk,lmj) /= jmr) cycle
      iwds = just(lri,lrk)
      iwdt = just(lrk,lrj)
      call prodab(3,jpel,iwds,iwdt,jwl,jwr,-wl,jper)
      ! st(2-3) ar(13)-c'(22)-bl(32)-
      if (jb_sys > 0) then
        iwds = just(lrk,lri)
        call prodab(3,jpel,iwds,iwdt,jwl,jwr,-wl1,jper)
      end if
    end do
  end do
end do

return

end subroutine st_head_dbl_tail_act

subroutine ts_head_dbl_tail_act(lra)

use gugaci_global, only: jb_sys, jml, jmr, jpel, jper, just, jwl, jwr, line, lrs, lsm_inn, norb_dz, norb_frz, w1, w1_ts
use Symmetry_Info, only: Mul
use Constants, only: Zero
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: lra
integer(kind=iwp) :: iwds, iwdt, lmas, lmi, lmij, lmj, lmk, lri, lrj, lrk, ni
real(kind=wp) :: vlop0, vlop1, w1ts1, w1ts2, w1ts3, w1ts4, wl

! ts(3-1) ar(23)-bl(20)-
! ts(3-2) (22)ar(23)-bl(31)-
! ts(3-2) ar(23)-c'(22)-bl(31)-
! ts(3-3) ar(23)-bl(31)-c"(22)-
! ts(3-4) ar(23)-bl(32)-c"(21)-
lmas = Mul(lsm_inn(lra),lsm_inn(lrs))
do lri=norb_frz+1,norb_dz
  lmi = lsm_inn(lri)
  do lrj=lri+1,norb_dz
    lmj = lsm_inn(lrj)
    lmij = Mul(lmi,lmj)
    if (lmas /= lmij) cycle
    w1ts1 = w1_ts(1)
    w1ts2 = w1_ts(2)
    w1ts3 = w1_ts(3)
    w1ts4 = w1_ts(4)
    ni = mod(lrj-lri,2)
    if (ni == 0) then
      w1ts1 = -w1ts1
      w1ts2 = -w1ts2
      w1ts3 = -w1ts3
      w1ts4 = -w1ts4
    end if
    !list = list3(lri,lrj,lrb)
    !-------------------------------------------------------------------
    ! ts(3-1) ar(23)-bl(20)-
    if ((jmr == 1) .and. (lmij == jml)) then
      iwdt = just(lri,lrj)
      iwds = just(lrj,lrj)
      vlop1 = w1*w1ts1             !d3-1
      vlop0 = Zero
      if (line == 26) then    !lri,lrj,lra
        call comp_loop(5,lri,lrj,lrs,lra,vlop0,vlop1,wl)
      end if
      if (line == 28) then    !lri,lrj,lrs,lra
        call comp_loop(7,lri,lrj,lrs,lra,vlop0,vlop1,wl)
      end if
      if (line == 29) then    !lri,lrj,lrs,lra
        call comp_loop(6,lri,lrj,lrs,lra,vlop0,vlop1,wl)
      end if
      !wl = -vlop1*(vint_ci(list))    !2.2 vlop0=0
      call prodab(3,jpel,iwdt,iwds,jwl,jwr,wl,jper)
    end if
    !-------------------------------------------------------------------
    vlop1 = w1*w1ts2             !d3-2
    vlop0 = Zero
    if (line == 26) then    !lri,lrj,lra
      call comp_loop(5,lri,lrj,lrs,lra,vlop0,vlop1,wl)
    end if
    if (line == 28) then    !lri,lrj,lrs,lra
      call comp_loop(7,lri,lrj,lrs,lra,vlop0,vlop1,wl)
    end if
    if (line == 29) then    !lri,lrj,lrs,lra
      call comp_loop(6,lri,lrj,lrs,lra,vlop0,vlop1,wl)
    end if
    !wl = -vlop1*(vint_ci(list))    !2.2 vlop0=0
    ! ts(3-2) (22)ar(23)-bl(31)-
    do lrk=norb_frz+1,lri-1
      lmk = lsm_inn(lrk)
      if ((Mul(lmk,lmi) /= jml) .or. (Mul(lmk,lmj) /= jmr)) cycle
      iwdt = just(lrk,lri)
      iwds = just(lrk,lrj)
      call prodab(3,jpel,iwdt,iwds,jwl,jwr,wl,jper)
    end do
    ! ts(3-2) ar(23)-c'(22)-bl(31)-
    do lrk=lri+1,lrj-1
      lmk = lsm_inn(lrk)
      if ((Mul(lmk,lmi) /= jml) .or. (Mul(lmk,lmj) /= jmr)) cycle
      iwdt = just(lri,lrk)
      iwds = just(lrk,lrj)
      call prodab(3,jpel,iwdt,iwds,jwl,jwr,-wl,jper)
    end do
    !-------------------------------------------------------------------
    ! ts(3-4) ar(23)-bl(32)-c"(21)-
    vlop1 = w1*w1ts4             !d3-4
    vlop0 = Zero
    if (line == 26) then    !lri,lrj,lra
      call comp_loop(5,lri,lrj,lrs,lra,vlop0,vlop1,wl)
    end if
    if (line == 28) then    !lri,lrj,lrs,lra
      call comp_loop(7,lri,lrj,lrs,lra,vlop0,vlop1,wl)
    end if
    if (line == 29) then    !lri,lrj,lrs,lra
      call comp_loop(6,lri,lrj,lrs,lra,vlop0,vlop1,wl)
    end if
    do lrk=lrj+1,norb_dz
      lmk = lsm_inn(lrk)
      if ((Mul(lmk,lmi) /= jml) .or. (Mul(lmk,lmj) /= jmr)) cycle
      iwdt = just(lri,lrk)
      iwds = just(lrj,lrk)
      call prodab(3,jpel,iwdt,iwds,jwl,jwr,wl,jper)
    end do
    !-------------------------------------------------------------------
    ! ts(3-3) ar(23)-bl(31)-c"(22)-
    if (jb_sys > 0) then
      vlop1 = w1*w1ts3             !d3-4
      vlop0 = Zero
      if (line == 26) then    !lri,lrj,lra
        call comp_loop(5,lri,lrj,lrs,lra,vlop0,vlop1,wl)
      end if
      if (line == 28) then    !lri,lrj,lrs,lra
        call comp_loop(7,lri,lrj,lrs,lra,vlop0,vlop1,wl)
      end if
      if (line == 29) then    !lri,lrj,lrs,lra
        call comp_loop(6,lri,lrj,lrs,lra,vlop0,vlop1,wl)
      end if
      !wl = -vlop1*vint_ci(list)    !2.2 vlop0=0
      do lrk=lrj+1,norb_dz
        lmk = lsm_inn(lrk)
        if ((Mul(lmk,lmi) /= jml) .or. (Mul(lmk,lmj) /= jmr)) cycle
        iwdt = just(lri,lrk)
        iwds = just(lrk,lrj)
        call prodab(3,jpel,iwdt,iwds,jwl,jwr,wl,jper)
      end do
    end if
  end do
end do

return

end subroutine ts_head_dbl_tail_act

subroutine stt_head_dbl_tail_act(lra)

use gugaci_global, only: jml, jmr, jpel, jper, just, jwl, jwr, line, lrs, lsm_inn, norb_dz, norb_frz, w1, w1_st1
use Symmetry_Info, only: Mul
use Constants, only: Zero
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: lra
integer(kind=iwp) :: iwdl, iwdr, lmi, lmij, lmj, lmk, lri, lrj, lrk, ni
real(kind=wp) :: vlop0, vlop1, w1st1, w1st2, w1st3, w1st4, wl

! st1(4-1) ar(01)-bl(31)-
! st1(4-2) (11)ar(23)-bl(31)-
! st1(4-3) ar(13)-c'(21)-bl(31)-
! st1(4-3) ar(13)-bl(31)-c"(21)-
! st1(4-4) ar(23)-c'(11)-bl(31)-
! st1(4-4) ar(23)-bl(31)-c"(11)-
do lri=norb_frz+1,norb_dz-1
  lmi = lsm_inn(lri)
  do lrj=lri+1,norb_dz
    lmj = lsm_inn(lrj)
    lmij = Mul(lmi,lmj)
    w1st1 = w1_st1(1)
    w1st2 = w1_st1(2)
    w1st3 = w1_st1(3)
    w1st4 = w1_st1(4)
    ni = mod(lrj-lri,2)
    if (ni == 0) then
      w1st1 = -w1st1
      w1st2 = -w1st2
      w1st3 = -w1st3
      w1st4 = -w1st4
    end if
    !list = list3(lri,lrj,lrb)
    if ((jml == 1) .and. (lmij == jmr)) then
      ! st1(4-1) ar(01)-bl(31)-
      iwdl = just(lri,lri)
      iwdr = just(lri,lrj)
      vlop1 = w1*w1st1
      vlop0 = Zero
      if (line == 26) then    !lri,lrj,lra
        call comp_loop(5,lri,lrj,lrs,lra,vlop0,vlop1,wl)
      end if
      if (line == 28) then    !lri,lrj,lrs,lra
        call comp_loop(7,lri,lrj,lrs,lra,vlop0,vlop1,wl)
      end if
      if (line == 29) then    !lri,lrj,lrs,lra
        call comp_loop(6,lri,lrj,lrs,lra,vlop0,vlop1,wl)
      end if
      !wl = -vlop1*(vint_ci(list))    !2.2 vlop0=0
      call prodab(3,jpel,iwdl,iwdr,jwl,jwr,wl,jper)
    end if
    ! st1(4-2) (11)ar(23)-bl(31)-
    vlop1 = w1*w1st2
    vlop0 = Zero
    if (line == 26) then    !lri,lrj,lra
      call comp_loop(5,lri,lrj,lrs,lra,vlop0,vlop1,wl)
    end if
    if (line == 28) then    !lri,lrj,lrs,lra
      call comp_loop(7,lri,lrj,lrs,lra,vlop0,vlop1,wl)
    end if
    if (line == 29) then    !lri,lrj,lrs,lra
      call comp_loop(6,lri,lrj,lrs,lra,vlop0,vlop1,wl)
    end if
    do lrk=norb_frz+1,lri-1
      lmk = lsm_inn(lrk)
      if ((Mul(lmk,lmi) /= jml) .or. (Mul(lmk,lmj) /= jmr)) cycle
      !wl = -vlop1*vint_ci(list)
      iwdl = just(lri,lrk)
      iwdr = just(lrk,lrj)
      call prodab(3,jpel,iwdl,iwdr,jwl,jwr,wl,jper)
    end do
    vlop1 = w1*w1st4
    vlop0 = Zero
    if (line == 26) then    !lri,lrj,lra
      call comp_loop(5,lri,lrj,lrs,lra,vlop0,vlop1,wl)
    end if
    if (line == 28) then    !lri,lrj,lrs,lra
      call comp_loop(7,lri,lrj,lrs,lra,vlop0,vlop1,wl)
    end if
    if (line == 29) then    !lri,lrj,lrs,lra
      call comp_loop(6,lri,lrj,lrs,lra,vlop0,vlop1,wl)
    end if
    ! st1(4-4) ar(23)-c'(11)-bl(31)-
    do lrk=lri+1,lrj-1
      lmk = lsm_inn(lrk)
      if ((Mul(lmi,lmk) /= jml) .or. (Mul(lmk,lmj) /= jmr)) cycle
      iwdl = just(lri,lrk)
      iwdr = just(lrk,lrj)
      call prodab(3,jpel,iwdl,iwdr,jwl,jwr,-wl,jper)
    end do
    ! st1(4-4) ar(23)-bl(31)-c"(11)-
    do lrk=lrj+1,norb_dz
      lmk = lsm_inn(lrk)
      if ((Mul(lmi,lmk) /= jml) .or. (Mul(lmj,lmk) /= jmr)) cycle
      iwdl = just(lri,lrk)
      iwdr = just(lrj,lrk)
      call prodab(3,jpel,iwdl,iwdr,jwl,jwr,wl,jper)
    end do

    vlop1 = w1*w1st3
    vlop0 = Zero
    if (line == 26) then    !lri,lrj,lra
      call comp_loop(5,lri,lrj,lrs,lra,vlop0,vlop1,wl)
    end if
    if (line == 28) then    !lri,lrj,lrs,lra
      call comp_loop(7,lri,lrj,lrs,lra,vlop0,vlop1,wl)
    end if
    if (line == 29) then    !lri,lrj,lrs,lra
      call comp_loop(6,lri,lrj,lrs,lra,vlop0,vlop1,wl)
    end if
    ! st1(4-3) ar(13)-c'(21)-bl(31)-
    do lrk=lri+1,lrj-1
      lmk = lsm_inn(lrk)
      if ((Mul(lmi,lmk) /= jml) .or. (Mul(lmk,lmj) /= jmr)) cycle
      iwdl = just(lrk,lri)
      iwdr = just(lrk,lrj)
      call prodab(3,jpel,iwdl,iwdr,jwl,jwr,-wl,jper)
    end do
    ! st1(4-3) ar(13)-bl(31)-c"(21)-
    do lrk=lrj+1,norb_dz
      lmk = lsm_inn(lrk)
      if ((Mul(lmi,lmk) /= jml) .or. (Mul(lmj,lmk) /= jmr)) cycle
      iwdl = just(lrk,lri)
      iwdr = just(lrj,lrk)
      call prodab(3,jpel,iwdl,iwdr,jwl,jwr,wl,jper)
    end do
  end do
end do

return

end subroutine stt_head_dbl_tail_act

subroutine tts_head_dbl_tail_act(lra)

use gugaci_global, only: jml, jmr, jpel, jper, just, jwl, jwr, line, lrs, lsm_inn, norb_dz, norb_frz, w1, w1_t1s
use Symmetry_Info, only: Mul
use Constants, only: Zero
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: lra
integer(kind=iwp) :: iwdl, iwdr, lmi, lmij, lmj, lmk, lri, lrj, lrk, ni
real(kind=wp) :: vlop0, vlop1, w1t1s1, w1t1s2, w1t1s3, w1t1s4, w1t1s5, w1t1s6, w1t1s7, wl

! t1s(5-1)   ar(13)-bl(10)-
! t1s(5-2)   ar(13)-bl(32)-
! t1s(5-2)   ar(13)-c'(11)-bl(32)-
! t1s(5-3)   ar(13)-bl(31)-c"(12)-
! t1s(5-4)   ar(13)-bl(32)-c"(11)-
! t1s(5-5)   drl(12)-
! t1s(5-6)   drl(12)-c"(12)-
! t1s(5-7)   drl(12)-c"(11)-
do lri=norb_frz+1,norb_dz
  lmi = lsm_inn(lri)
  do lrj=lri+1,norb_dz
    lmj = lsm_inn(lrj)
    lmij = Mul(lmi,lmj)
    w1t1s1 = w1_t1s(1)
    w1t1s2 = w1_t1s(2)
    w1t1s3 = w1_t1s(3)
    w1t1s4 = w1_t1s(4)
    w1t1s5 = w1_t1s(5)
    w1t1s6 = w1_t1s(6)
    w1t1s7 = w1_t1s(7)
    ni = mod(lrj-lri,2)
    if (ni == 0) then
      w1t1s1 = -w1t1s1
      w1t1s2 = -w1t1s2
      w1t1s3 = -w1t1s3
      w1t1s4 = -w1t1s4
      !w1t1s5 = -w1t1s5
      !w1t1s6 = -w1t1s6
      !w1t1s7 = -w1t1s7
    end if
    !list = list3(lri,lrj,lrb)
    if ((jml == lmij) .and. (jmr == 1)) then
      ! t1s(5-1)   ar(13)-bl(10)-
      vlop1 = w1*w1t1s1
      vlop0 = Zero
      if (line == 26) then    !lri,lrj,lra
        call comp_loop(5,lri,lrj,lrs,lra,vlop0,vlop1,wl)
      end if
      if (line == 28) then    !lri,lrj,lrs,lra
        call comp_loop(7,lri,lrj,lrs,lra,vlop0,vlop1,wl)
      end if
      if (line == 29) then    !lri,lrj,lrs,lra
        call comp_loop(6,lri,lrj,lrs,lra,vlop0,vlop1,wl)
      end if
      !wl = -vlop1*vint_ci(list)
      iwdl = just(lri,lrj)
      iwdr = just(lrj,lrj)
      call prodab(3,jpel,iwdl,iwdr,jwl,jwr,wl,jper)
    end if
    vlop1 = w1*w1t1s2
    vlop0 = Zero
    if (line == 26) then    !lri,lrj,lra
      call comp_loop(5,lri,lrj,lrs,lra,vlop0,vlop1,wl)
    end if
    if (line == 28) then    !lri,lrj,lrs,lra
      call comp_loop(7,lri,lrj,lrs,lra,vlop0,vlop1,wl)
    end if
    if (line == 29) then    !lri,lrj,lrs,lra
      call comp_loop(6,lri,lrj,lrs,lra,vlop0,vlop1,wl)
    end if
    !wl = -vlop1*vint_ci(list)
    do lrk=norb_frz+1,lri-1
      lmk = lsm_inn(lrk)
      ! t1s(5-2)   (11)ar(13)-bl(32)-
      if ((jml /= Mul(lmk,lmi)) .or. (jmr /= Mul(lmk,lmj))) cycle
      !vlop1 = w1*w1t1s2
      !wl = -vlop1*vint_ci(list)
      iwdl = just(lrk,lri)
      iwdr = just(lrj,lrk)
      call prodab(3,jpel,iwdl,iwdr,jwl,jwr,wl,jper)
    end do
    do lrk=lri+1,lrj-1
      lmk = lsm_inn(lrk)
      ! t1s(5-2)   ar(13)-c'(11)-bl(32)-
      if ((jml /= Mul(lmi,lmk)) .or. (jmr /= Mul(lmk,lmj))) cycle
      !vlop1 = w1*w1t1s2
      !wl = -vlop1*vint_ci(list)
      iwdl = just(lri,lrk)
      iwdr = just(lrj,lrk)
      call prodab(3,jpel,iwdl,iwdr,jwl,jwr,-wl,jper)
    end do
    vlop1 = w1*w1t1s3
    vlop0 = Zero
    if (line == 26) then    !lri,lrj,lra
      call comp_loop(5,lri,lrj,lrs,lra,vlop0,vlop1,wl)
    end if
    if (line == 28) then    !lri,lrj,lrs,lra
      call comp_loop(7,lri,lrj,lrs,lra,vlop0,vlop1,wl)
    end if
    if (line == 29) then    !lri,lrj,lrs,lra
      call comp_loop(6,lri,lrj,lrs,lra,vlop0,vlop1,wl)
    end if
    do lrk=lrj+1,norb_dz
      lmk = lsm_inn(lrk)
      ! t1s(5-3)   ar(13)-bl(31)-c"(12)-
      if ((jml /= Mul(lmi,lmk)) .or. (jmr /= Mul(lmj,lmk))) cycle
      !vlop1 = w1*w1t1s3
      !wl = -vlop1*vint_ci(list)
      iwdl = just(lri,lrk)
      iwdr = just(lrk,lrj)
      call prodab(3,jpel,iwdl,iwdr,jwl,jwr,wl,jper)
    end do
    vlop1 = w1*w1t1s4
    vlop0 = Zero
    if (line == 26) then    !lri,lrj,lra
      call comp_loop(5,lri,lrj,lrs,lra,vlop0,vlop1,wl)
    end if
    if (line == 28) then    !lri,lrk,lrs,lra
      call comp_loop(7,lri,lrj,lrs,lra,vlop0,vlop1,wl)
    end if
    if (line == 29) then    !lri,lrk,lrs,lra
      call comp_loop(6,lri,lrj,lrs,lra,vlop0,vlop1,wl)
    end if
    !wl = -vlop1*vint_ci(list)
    do lrk=lrj+1,norb_dz
      lmk = lsm_inn(lrk)
      ! t1s(5-4)   ar(13)-bl(32)-c"(11)-
      if ((jml /= Mul(lmi,lmk)) .or. (jmr /= Mul(lmj,lmk))) cycle
      iwdl = just(lri,lrk)
      iwdr = just(lrj,lrk)
      call prodab(3,jpel,iwdl,iwdr,jwl,jwr,wl,jper)
    end do
    if ((lmij == jml) .and. (lmij == jmr)) then
      ! t1s(5-5)   (11)drl(12)-
      vlop1 = w1*w1t1s5
      vlop0 = Zero
      if (line == 26) then    !lrj,lra
        call comp_loop(9,lrj,0,0,lra,vlop0,vlop1,wl)
      end if
      if (line == 28) then    !lrj,lrs,lra
        call comp_loop(12,lrj,0,lrs,lra,vlop0,vlop1,wl)
      end if
      if (line == 29) then    !lrj,lrs,lra
        call comp_loop(11,lrj,0,lrs,lra,vlop0,vlop1,wl)
      end if
      !wl = -vlop1*voint(lrb,lrj)
      iwdl = just(lri,lrj)
      iwdr = just(lrj,lri)
      call prodab(3,jpel,iwdl,iwdr,jwl,jwr,wl,jper)
      ! t1s(5-6)   drl(11)-c"(12)-
      vlop1 = w1*w1t1s6
      vlop0 = Zero
      if (line == 26) then    !lri,lra
        call comp_loop(9,lri,0,0,lra,vlop0,vlop1,wl)
      end if
      if (line == 28) then    !lri,lrs,lra
        call comp_loop(12,lri,0,lrs,lra,vlop0,vlop1,wl)
      end if
      if (line == 29) then    !lri,lrs,lra
        call comp_loop(11,lri,0,lrs,lra,vlop0,vlop1,wl)
      end if
      !wl = -vlop1*voint(lrb,lri)
      iwdl = just(lri,lrj)
      iwdr = just(lrj,lri)
      call prodab(3,jpel,iwdl,iwdr,jwl,jwr,wl,jper)
      ! t1s(5-7)   drl(12)-c"(11)-
      vlop1 = w1*w1t1s7
      vlop0 = Zero
      if (line == 26) then    !lri,lra
        call comp_loop(9,lri,0,0,lra,vlop0,vlop1,wl)
      end if
      if (line == 28) then    !lri,lrs,lra
        call comp_loop(12,lri,0,lrs,lra,vlop0,vlop1,wl)
      end if
      if (line == 29) then    !lri,lrs,lra
        call comp_loop(11,lri,0,lrs,lra,vlop0,vlop1,wl)
      end if
      !wl = -vlop1*voint(lrb,lri)
      iwdl = just(lri,lrj)
      iwdr = iwdl
      call prodab(3,jpel,iwdl,iwdr,jwl,jwr,wl,jper)
    end if
  end do
end do

return

end subroutine tts_head_dbl_tail_act

subroutine tt_head_dbl_tail_act(lra)

use gugaci_global, only: jml, jmr, jpad, jpadl, jpel, jper, just, jwl, jwr, line, lrs, lsm_inn, norb_dz, norb_frz, w0, w0_tt, w1, &
                         w1_tt
use Symmetry_Info, only: Mul
use Constants, only: Zero
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: lra
integer(kind=iwp) :: iwdl, iwdr, jmlr, lmi, lmij, lmj, lmk, lmki, lmkj, lri, lrj, lrk, ni
real(kind=wp) :: vlop0, vlop1, w0tt1, w0tt2, w0tt3, w1tt1, w1tt2, wl, wltmp

! tt(11-1) (22)ar(23)-bl(32)-
! tt(11-1) ar(23)-c'(22)-bl(32)-
! tt(11-1) ar(23)-bl(32)-c"(22)-
! tt(11-2) (22)drl(22)-
! tt(11-2) drl(22)-c"(22)-
! tt(11-3) (22)drl(33)-
! tt(11-3) drl(33)-c"(22)-
! tt(11-3) drl(33)-c"(22)-c"(22)-
do lri=norb_frz+1,norb_dz-1
  lmi = lsm_inn(lri)
  do lrj=lri+1,norb_dz
    lmj = lsm_inn(lrj)
    lmij = Mul(lmi,lmj)
    jmlr = Mul(jml,jmr)
    if (lmij /= jmlr) cycle
    w0tt1 = w0_tt(1)
    w1tt1 = w1_tt(1)
    ni = mod(lrj-lri,2)
    if (ni == 0) then
      w0tt1 = -w0tt1
      w1tt1 = -w1tt1
    end if
    !list = list3(lri,lrj,lrb)
    vlop0 = w0*w0tt1
    vlop1 = w1*w1tt1
    if (line == 26) then    !lri,lrj,lra
      call comp_loop(5,lri,lrj,lrs,lra,vlop0,vlop1,wl)
    end if
    if (line == 28) then    !lri,lrj,lrs,lra
      call comp_loop(7,lri,lrj,lrs,lra,vlop0,vlop1,wl)
    end if
    if (line == 29) then    !lri,lrj,lrs,lra
      call comp_loop(6,lri,lrj,lrs,lra,vlop0,vlop1,wl)
    end if
    !wl = (vlop0-vlop1)*vint_ci(list)-2*vlop0*vint_ci(list+1)
    ! tt(11-1) (22)ar(23)-bl(32)-
    do lrk=norb_frz+1,lri-1
      lmk = lsm_inn(lrk)
      lmki = Mul(lmk,lmi)
      lmkj = Mul(lmk,lmj)
      if ((lmki == jml) .and. (lmkj == jmr)) then
        iwdl = just(lrk,lri)
        iwdr = just(lrk,lrj)
        call prodab(3,jpel,iwdl,iwdr,jwl,jwr,wl,jper)
      end if
    end do
    ! tt(11-1) ar(23)-bl(32)-c"(22)-    act -c"-
    do lrk=lrj+1,norb_dz
      lmk = lsm_inn(lrk)
      lmki = Mul(lmk,lmi)
      lmkj = Mul(lmk,lmj)
      if ((lmki == jml) .and. (lmkj == jmr)) then
        iwdl = just(lri,lrk)
        iwdr = just(lrj,lrk)
        call prodab(3,jpel,iwdl,iwdr,jwl,jwr,wl,jper)
      end if
    end do
    ! tt(11-1) ar(23)-c'(22)-bl(32)-    act -c"-
    do lrk=lri+1,lrj-1
      lmk = lsm_inn(lrk)
      lmki = Mul(lmk,lmi)
      lmkj = Mul(lmk,lmj)
      if ((lmki == jml) .and. (lmkj == jmr)) then
        iwdl = just(lri,lrk)
        iwdr = just(lrk,lrj)
        call prodab(3,jpel,iwdl,iwdr,jwl,jwr,-wl,jper)
      end if
    end do
  end do
end do

if (jpad /= jpadl) return
if (jwl >= jwr) return

w0tt2 = w0_tt(2)
w1tt2 = w1_tt(2)
w0tt3 = w0_tt(3)

do lri=norb_frz+1,norb_dz-1
  lmi = lsm_inn(lri)
  do lrj=lri+1,norb_dz
    lmj = lsm_inn(lrj)
    lmij = Mul(lmi,lmj)
    if ((lmij /= jml) .or. (lmij /= jmr)) cycle
    if (jwl >= jwr) cycle
    ! tt(11-2) (22)drl(22)-
    ! tt(11-2) drl(22)-c"(22)-
    iwdl = just(lri,lrj)
    iwdr = iwdl
    vlop0 = w0*w0tt2
    vlop1 = w1*w1tt2
    if (line == 26) then    !lrj,lra
      call comp_loop(9,lrj,0,0,lra,vlop0,vlop1,wl)
    end if
    if (line == 28) then    !lrj,lrs,lra
      call comp_loop(12,lrj,0,lrs,lra,vlop0,vlop1,wl)
    end if
    if (line == 29) then    !lrj,lrs,lra
      call comp_loop(11,lrj,0,lrs,lra,vlop0,vlop1,wl)
    end if
    if (line == 26) then    !lri,lra
      call comp_loop(9,lri,0,0,lra,vlop0,vlop1,wltmp)
    end if
    if (line == 28) then    !lri,lrs,lra
      call comp_loop(12,lri,0,lrs,lra,vlop0,vlop1,wltmp)
    end if
    if (line == 29) then    !lri,lrs,lra
      call comp_loop(11,lri,0,lrs,lra,vlop0,vlop1,wltmp)
    end if
    wl = wl+wltmp
    !wl = (vlop0-vlop1)*(voint(lrb,lri)+voint(lrb,lrj))
    vlop0 = w0*w0tt3
    vlop1 = Zero
    do lrk=1,norb_dz
      if (lrk == lri) cycle
      if (lrk == lrj) cycle
      ! tt(11-3) drl(33)-c"(22)-c"(22)-
      ! tt(11-3) (22)drl(33)-c"(22)-
      ! tt(11-3) (22)(22)drl(33)-
      if (line == 26) then    !lrk,lra
        call comp_loop(9,lrk,0,0,lra,vlop0,vlop1,wltmp)
      end if
      if (line == 28) then    !lrk,lrs,lra
        call comp_loop(12,lrk,0,lrs,lra,vlop0,vlop1,wltmp)
      end if
      if (line == 29) then    !lrk,lrs,lra
        call comp_loop(11,lrk,0,lrs,lra,vlop0,vlop1,wltmp)
      end if
      wl = wl+wltmp
    end do
    call prodab(3,jpel,iwdl,iwdr,jwl,jwr,wl,jper)
  end do
end do

return

end subroutine tt_head_dbl_tail_act

subroutine tttt_head_dbl_tail_act(lra)

use gugaci_global, only: jml, jmr, jpad, jpadl, jpel, jper, just, jwl, jwr, line, lrs, lsm_inn, norb_dz, norb_frz, w0, w0_t1t1, &
                         w1, w1_t1t1
use Symmetry_Info, only: Mul
use Constants, only: Zero
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: lra
integer(kind=iwp) :: iwdl, iwdr, lmi, lmij, lmj, lmm, lri, lrj, lrk, lrm, ni
real(kind=wp) :: vlop0, vlop1, w0tt1, w1tt1, wl, wltmp

! t1t1(12-1)  ar(13)-bl(31)-
! t1t1(12-1)  ar(13)-c'(11)-bl(31)-
! t1t1(12-1)  ar(13)-bl(31)-c"(11)-
do lri=norb_frz+1,norb_dz-1
  lmi = lsm_inn(lri)
  do lrj=lri+1,norb_dz
    lmj = lsm_inn(lrj)
    w0tt1 = w0_t1t1(1)
    w1tt1 = w1_t1t1(1)
    ni = mod(lrj-lri,2)
    if (ni == 0) then
      w0tt1 = -w0tt1
      w1tt1 = -w1tt1
    end if
    !list = list4(lri,lrj,lra,lrb)
    vlop0 = w0*w0tt1
    vlop1 = w1*w1tt1
    if (line == 26) then    !lri,lrj,lra
      call comp_loop(5,lri,lrj,lrs,lra,vlop0,vlop1,wl)
    end if
    if (line == 28) then    !lri,lrj,lrs,lra
      call comp_loop(7,lri,lrj,lrs,lra,vlop0,vlop1,wl)
    end if
    if (line == 29) then    !lri,lrj,lrs,lra
      call comp_loop(6,lri,lrj,lrs,lra,vlop0,vlop1,wl)
    end if
    !wl = (vlop0-vlop1)*vint_ci(list)-2*vlop0*vint_ci(list+1)
    ! t1t1(12-1)  (11)ar(13)-bl(31)-
    do lrm=norb_frz+1,lri-1
      lmm = lsm_inn(lrm)
      if ((jml /= Mul(lmi,lmm)) .or. (jmr /= Mul(lmm,lmj))) cycle
      iwdl = just(lrm,lri)
      iwdr = just(lrm,lrj)
      call prodab(3,jpel,iwdl,iwdr,jwl,jwr,wl,jper)
    end do
    ! t1t1(12-1)  ar(13)-bl(31)-c"(11)-
    do lrm=lrj+1,norb_dz
      lmm = lsm_inn(lrm)
      if ((jml /= Mul(lmi,lmm)) .or. (jmr /= Mul(lmm,lmj))) cycle
      iwdl = just(lri,lrm)
      iwdr = just(lrj,lrm)
      call prodab(3,jpel,iwdl,iwdr,jwl,jwr,wl,jper)
    end do
    ! t1t1(12-1)  ar(13)-c'(11)-bl(31)-
    do lrm=lri+1,lrj-1
      lmm = lsm_inn(lrm)
      if ((jml /= Mul(lmi,lmm)) .or. (jmr /= Mul(lmm,lmj))) cycle
      iwdl = just(lri,lrm)
      iwdr = just(lrm,lrj)
      call prodab(3,jpel,iwdl,iwdr,jwl,jwr,-wl,jper)
    end do
  end do
end do

if (jpad /= jpadl) return
if (jwl >= jwr) return

do lri=norb_frz+1,norb_dz
  lmi = lsm_inn(lri)
  do lrj=lri+1,norb_dz               !bbs_tmp
    lmj = lsm_inn(lrj)
    lmij = Mul(lmi,lmj)
    if (jml /= lmij) cycle
    vlop0 = w0*w0_t1t1(2)
    vlop1 = w1*w1_t1t1(2)
    ! t1t1(12-2)  (11)drl(11)-
    if (line == 26) then    !lri,lra
      call comp_loop(9,lri,0,0,lra,vlop0,vlop1,wl)
    end if
    if (line == 28) then    !lrk,lrs,lra
      call comp_loop(12,lri,0,lrs,lra,vlop0,vlop1,wl)
    end if
    if (line == 29) then    !lrk,lrs,lra
      call comp_loop(11,lri,0,lrs,lra,vlop0,vlop1,wl)
    end if
    ! t1t1(12-2)  drl(11)-c"(11)-
    if (line == 26) then    !lrj,lra
      call comp_loop(9,lrj,0,0,lra,vlop0,vlop1,wltmp)
    end if
    if (line == 28) then    !lrj,lrs,lra
      call comp_loop(12,lrj,0,lrs,lra,vlop0,vlop1,wltmp)
    end if
    if (line == 29) then    !lrj,lrs,lra
      call comp_loop(11,lrj,0,lrs,lra,vlop0,vlop1,wltmp)
    end if
    wl = wl+wltmp
    ! t1t1(12-3)  (11)(11)drl(33)-
    ! t1t1(12-3)  (11)drl(33)-c"(11)-
    ! t1t1(12-3)  drl(33)-c"(11)-c"(11)-
    do lrk=1,norb_dz
      if (lrk == lri) cycle
      if (lrk == lrj) cycle
      vlop0 = w0*w0_t1t1(3)
      vlop1 = Zero
      if (line == 26) then    !lrk,lra
        call comp_loop(9,lrk,0,0,lra,vlop0,vlop1,wltmp)
      end if
      if (line == 28) then    !lrk,lrs,lra
        call comp_loop(12,lrk,0,lrs,lra,vlop0,vlop1,wltmp)
      end if
      if (line == 29) then    !lrk,lrs,lra
        call comp_loop(11,lrk,0,lrs,lra,vlop0,vlop1,wltmp)
      end if
      wl = wl+wltmp
    end do
    iwdl = just(lri,lrj)
    iwdr = iwdl
    call prodab(3,jpel,iwdl,iwdr,jwl,jwr,wl,jper)
  end do
end do

return

end subroutine tttt_head_dbl_tail_act

subroutine dd_head_dbl_tail_act(lra)

use gugaci_global, only: jml, jmr, jpel, jper, jud, jwl, jwr, line, lrs, lsm_inn, norb_dz, norb_frz, w0, w0_dd, w1, w1_dd
use Constants, only: Zero
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: lra
integer(kind=iwp) :: imil, imir, iwdl, iwdr, lril, lrir, lrk, ni
real(kind=wp) :: vlop0, vlop1, w0dd1, w1dd1, wl, wltmp

! dd(19-1) ar(23)-bl(32)-
! dd(19-2) drl(22)-
! dd(19-3) drl(33)-
! dd(19-3) drl(33)-c"(22)-
do lril=norb_frz+1,norb_dz
  imil = lsm_inn(lril)
  if (imil /= jml) cycle
  iwdl = jud(lril)
  do lrir=lril,norb_dz
    imir = lsm_inn(lrir)
    if (imir /= jmr) cycle
    iwdr = jud(lrir)

    w0dd1 = w0_dd(1)
    w1dd1 = w1_dd(1)
    ni = mod(lrir-lril,2)
    if (ni == 0) then
      w0dd1 = -w0dd1
      w1dd1 = -w1dd1
    end if

    if ((lril == lrir) .and. (jwl < jwr)) then
      ! dd(19-2) drl(22)-
      vlop0 = w0*w0_dd(2)
      vlop1 = w1*w1_dd(2)
      if (line == 26) then    !lril,lra
        call comp_loop(9,lril,0,0,lra,vlop0,vlop1,wl)
      end if
      if (line == 28) then    !lril,lrs,lra
        call comp_loop(12,lril,0,lrs,lra,vlop0,vlop1,wl)
      end if
      if (line == 29) then    !lril,lrs,lra
        call comp_loop(11,lril,0,lrs,lra,vlop0,vlop1,wl)
      end if
      !wl = (vlop0-vlop1)*voint(lra,lril)
      !call prodab(3,jpel,iwdl,iwdr,jwl,jwr,wl,jper)
      ! dd(19-3) drl(33)-
      vlop0 = w0*w0_dd(3)
      vlop1 = Zero
      do lrk=1,norb_dz
        if (lrk == lril) cycle
        !list = list3(lrs,lra,lrk)
        if (line == 26) then    !lrk,lra
          call comp_loop(9,lrk,0,0,lra,vlop0,vlop1,wltmp)     !wyb
        end if
        if (line == 28) then    !lrk,lrs,lra
          call comp_loop(12,lrk,0,lrs,lra,vlop0,vlop1,wltmp)  !wyb
        end if
        if (line == 29) then    !lrk,lrs,lra
          call comp_loop(11,lrk,0,lrs,lra,vlop0,vlop1,wltmp)  !wyb
        end if
        !wl = wl+vlop0*voint(lrk,lra)
        wl = wltmp+wl
      end do
      call prodab(3,jpel,iwdl,iwdr,jwl,jwr,wl,jper)
    end if
    if (lril /= lrir) then
      ! dd(19-1) ar(23)-bl(32)-
      vlop0 = w0*w0dd1
      vlop1 = w1*w1dd1
      if (line == 26) then   !lril,lrir,lra
        call comp_loop(5,lril,lrir,lrs,lra,vlop0,vlop1,wl)
      end if
      if (line == 28) then   !lril,lrir,lrs,lra
        call comp_loop(7,lril,lrir,lrs,lra,vlop0,vlop1,wl)
      end if
      if (line == 29) then   !lril,lrir,lrs,lra
        call comp_loop(6,lril,lrir,lrs,lra,vlop0,vlop1,wl)
      end if
      !list = list4(lril,lrir,lrs,lra)
      !wl = vlop0*(vint_ci(list)-2*vint_ci(list+1))-vlop1*vint_ci(list) !1.1       !!!!!
      call prodab(3,jpel,iwdl,iwdr,jwl,jwr,wl,jper)
    end if
  end do
end do

return

end subroutine dd_head_dbl_tail_act

subroutine dddd_head_dbl_tail_act(lra)

use gugaci_global, only: jml, jmr, jpel, jper, jud, jwl, jwr, line, lrs, lsm_inn, norb_dz, norb_frz, w0, w0_d1d1, w1, w1_d1d1
use Constants, only: Zero
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: lra
integer(kind=iwp) :: imil, imir, iwdl, iwdr, lril, lrir, lrk, ni
real(kind=wp) :: vlop0, vlop1, w0dd1, w1dd1, wl, wltmp

! d1d1(20-1) ar(13)-bl(31)-
! d1d1(20-1) drl(11)-
! d1d1(20-1) drl(33)-
! d1d1(20-1) drl(33)-c"(11)-
do lril=norb_frz+1,norb_dz
  imil = lsm_inn(lril)
  if (imil /= jml) cycle
  iwdl = jud(lril)
  do lrir=lril,norb_dz
    imir = lsm_inn(lrir)
    if (imir /= jmr) cycle
    w0dd1 = w0_d1d1(1)
    w1dd1 = w1_d1d1(1)
    ni = mod(lrir-lril,2)
    if (ni == 0) then
      w0dd1 = -w0dd1
      w1dd1 = -w1dd1
    end if
    iwdr = jud(lrir)
    if ((lril == lrir) .and. (jwl < jwr)) then
      ! d1d1(20-1) drl(11)-
      vlop0 = w0*w0_d1d1(2)
      vlop1 = w1*w1_d1d1(2)
      if (line == 26) then    !lril,lra
        call comp_loop(9,lril,0,0,lra,vlop0,vlop1,wl)
      end if
      if (line == 28) then    !lril,lrs,lra
        call comp_loop(12,lril,0,lrs,lra,vlop0,vlop1,wl)
      end if
      if (line == 29) then    !lril,lrs,lra
        call comp_loop(11,lril,0,lrs,lra,vlop0,vlop1,wl)
      end if
      ! d1d1(20-1) drl(33)-
      ! d1d1(20-1) drl(33)-c"(11)-
      vlop0 = w0*w0_d1d1(3)
      vlop1 = Zero
      do lrk=1,norb_dz
        if (lrk == lril) cycle
        if (line == 26) then    !lrk,lra
          call comp_loop(9,lrk,0,0,lra,vlop0,vlop1,wltmp)
        end if
        if (line == 28) then    !lrk,lrs,lra
          call comp_loop(12,lrk,0,lrs,lra,vlop0,vlop1,wltmp)
        end if
        if (line == 29) then    !lrk,lrs,lra
          call comp_loop(11,lrk,0,lrs,lra,vlop0,vlop1,wltmp)
        end if
        wl = wl+wltmp
      end do
      call prodab(3,jpel,iwdl,iwdr,jwl,jwr,wl,jper)
    end if
    if (lril /= lrir) then
      ! d1d1(20-1) ar(13)-bl(31)-
      vlop0 = w0*w0dd1
      vlop1 = w1*w1dd1
      if (line == 26) then   !lril,lrir,lra
        call comp_loop(5,lril,lrir,lrs,lra,vlop0,vlop1,wl)
      end if
      if (line == 28) then   !lril,lrir,lrs,lra
        call comp_loop(7,lril,lrir,lrs,lra,vlop0,vlop1,wl)
      end if
      if (line == 29) then   !lril,lrir,lrs,lra
        call comp_loop(6,lril,lrir,lrs,lra,vlop0,vlop1,wl)
      end if
      !list = list3(lril,lrir,lra)
      !wl = (vlop0-vlop1)*vint_ci(list)-2*vlop0*vint_ci(list+1)   !2
      call prodab(3,jpel,iwdl,iwdr,jwl,jwr,wl,jper)
    end if
  end do
end do

return

end subroutine dddd_head_dbl_tail_act

subroutine dd1_head_dbl_tail_act(lra)

use gugaci_global, only: jml, jmr, jpel, jper, jud, jwl, jwr, line, lrs, lsm_inn, norb_dz, norb_frz, w1, w1_dd1
use Constants, only: Zero
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: lra
integer(kind=iwp) :: iwdl, iwdr, lmi, lmj, lri, lrj, ni
real(kind=wp) :: vlop0, vlop1, wl

do lri=norb_frz+1,norb_dz
  lmi = lsm_inn(lri)
  iwdl = jud(lri)
  do lrj=lri+1,norb_dz
    ! dd1(21) ar(23)-bl(31)-
    lmj = lsm_inn(lrj)
    if ((jml /= lmi) .or. (jmr /= lmj)) cycle
    vlop1 = w1*w1_dd1
    vlop0 = Zero
    ni = mod(lrj-lri,2)
    if (ni == 0) then
      vlop1 = -vlop1
    end if

    if (line == 26) then   !lri,lrj,lra
      call comp_loop(5,lri,lrj,lrs,lra,vlop0,vlop1,wl)
    end if
    if (line == 28) then   !lril,lrj,lrs,lra
      call comp_loop(7,lri,lrj,lrs,lra,vlop0,vlop1,wl)
    end if
    if (line == 29) then   !lri,lrj,lrs,lra
      call comp_loop(6,lri,lrj,lrs,lra,vlop0,vlop1,wl)
    end if
    !list = list3(lri,lrj,lra)
    !wl = -vlop1*vint_ci(list)
    iwdr = jud(lrj)
    call prodab(3,jpel,iwdl,iwdr,jwl,jwr,wl,jper)
  end do
end do

return

end subroutine dd1_head_dbl_tail_act

subroutine d1d_head_dbl_tail_act(lra)

use gugaci_global, only: jml, jmr, jpel, jper, jud, jwl, jwr, line, lrs, lsm_inn, norb_dz, norb_frz, w1, w1_d1d
use Constants, only: Zero
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: lra
integer(kind=iwp) :: iwdl, iwdr, lmi, lmj, lri, lrj
real(kind=wp) :: vlop0, vlop1, wl

! d1d(22-1)   ar(13)-bl(32)-
! d1d(22-2)   drl(12)-
do lri=norb_frz+1,norb_dz-1
  lmi = lsm_inn(lri)
  iwdl = jud(lri)
  do lrj=lri+1,norb_dz
    lmj = lsm_inn(lrj)
    if ((jml /= lmi) .or. (jmr /= lmj)) cycle
    vlop1 = w1*w1_d1d(1)
    vlop0 = Zero
    if (mod(lrj-lri,2) == 0) then
      vlop1 = -vlop1
    end if
    if (line == 26) then   !lri,lrj,lra
      call comp_loop(5,lri,lrj,lrs,lra,vlop0,vlop1,wl)
    end if
    if (line == 28) then   !lril,lrj,lrs,lra
      call comp_loop(7,lri,lrj,lrs,lra,vlop0,vlop1,wl)
    end if
    if (line == 29) then   !lri,lrj,lrs,lra
      call comp_loop(6,lri,lrj,lrs,lra,vlop0,vlop1,wl)
    end if
    !list = list3(lri,lrj,lra)
    !wl = -vlop1*vint_ci(list)
    iwdr = jud(lrj)
    call prodab(3,jpel,iwdl,iwdr,jwl,jwr,wl,jper)
  end do
end do
vlop1 = w1*w1_d1d(2)
vlop0 = Zero
do lri=norb_frz+1,norb_dz
  ! d1d(22-2)   drl(12)-
  lmi = lsm_inn(lri)
  if ((jml /= lmi) .or. (jmr /= lmi)) cycle
  if (line == 26) then    !lri,lra
    call comp_loop(9,lri,0,0,lra,vlop0,vlop1,wl)
  end if
  if (line == 28) then    !lri,lrs,lra
    call comp_loop(12,lri,0,lrs,lra,vlop0,vlop1,wl)
  end if
  if (line == 29) then    !lri,lrs,lra
    call comp_loop(11,lri,0,lrs,lra,vlop0,vlop1,wl)
  end if
  !vlop1 = w1*w1_d1d(2)
  !wl = -vlop1*voint(lri,lra)
  iwdl = jud(lri)
  iwdr = iwdl
  call prodab(3,jpel,iwdl,iwdr,jwl,jwr,wl,jper)
end do

return

end subroutine d1d_head_dbl_tail_act

subroutine sv_head_dbl_tail_act(lra)

use gugaci_global, only: jb_sys, jml, jpel, jper, just, jwl, jwr, line, lrs, lsm_inn, norb_dz, norb_frz, w0, w0_sv, w1, w1_sv
use Symmetry_Info, only: Mul
use Constants, only: Zero
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: lra
integer(kind=iwp) :: iwdl, iwdr, lmi, lmij, lmj, lri, lrj, ni
real(kind=wp) :: vlop0, vlop1, w0sv1, w0sv2, w1sv1, w1sv2, wl

! sv(10-1) ar(13)-br(23)-
! sv(10-2) ar(23)-br(13)-
! sv(10-3) drr(03)-
iwdr = 0
do lri=norb_frz+1,norb_dz
  lmi = lsm_inn(lri)
  do lrj=lri,norb_dz
    lmj = lsm_inn(lrj)
    lmij = Mul(lmi,lmj)
    if (lmij /= jml) cycle
    w0sv1 = w0_sv(1)
    w1sv1 = w1_sv(1)
    w0sv2 = w0_sv(2)
    w1sv2 = w1_sv(2)
    ni = mod(lrj-lri,2)
    if (ni == 0) then
      w0sv1 = -w0sv1
      w1sv1 = -w1sv1
      w0sv2 = -w0sv2
      w1sv2 = -w1sv2
    end if
    iwdl = just(lri,lrj)
    if (lri == lrj) then
      vlop0 = w0*w0_sv(3)            !d10-3
      vlop1 = Zero
      !wl = vlop0*voint(lra,lri)/2
      if (line == 25) then    !lri,lra
        ! wyb call comp_loop(8,lri,lrg,lrs,lra,vlop0,vlop1,wl)
        call comp_loop(8,lri,0,0,lra,vlop0,vlop1,wl)
      end if
      if (line == 27) then    !lri,lrs,lra
        ! wyb call comp_loop(10,lri,lrs,lra,lra,vlop0,vlop1,wl)
        call comp_loop(10,lri,0,lrs,lra,vlop0,vlop1,wl)
      end if
      call prodab(3,jpel,iwdl,iwdr,jwl,jwr,wl,jper)
    else
      vlop0 = w0*w0sv2             !d10-2
      vlop1 = w1*w1sv2
      !list = list3(lri,lrj,lra)
      !wl = (vlop0+vlop1)*vint_ci(list)        !2.1          !!!!!
      if (line == 25) then    !lri,lrj,lra
        call comp_loop(3,lri,lrj,lrs,lra,vlop0,vlop1,wl)
      end if
      if (line == 27) then    !lri,lrj,lrs,lra
        call comp_loop(4,lri,lrj,lrs,lra,vlop0,vlop1,wl)
      end if
      call prodab(3,jpel,iwdl,iwdr,jwl,jwr,wl,jper)
      if (jb_sys > 0) then
        iwdl = just(lrj,lri)
        vlop0 = w0*w0sv1
        vlop1 = w1*w1sv1
        !list = list3(lri,lrj,lra)
        !wl = (vlop0+vlop1)*vint_ci(list)        !2.1          !!!!!
        if (line == 25) then    !lri,lrj,lra
          call comp_loop(3,lri,lrj,lrs,lra,vlop0,vlop1,wl)
        end if
        if (line == 27) then    !lri,lrj,lrs,lra
          call comp_loop(4,lri,lrj,lrs,lra,vlop0,vlop1,wl)
        end if
        call prodab(3,jpel,iwdl,iwdr,jwl,jwr,wl,jper)
      end if
    end if
  end do
end do

return

end subroutine sv_head_dbl_tail_act

subroutine sd_head_dbl_tail_act(lra,lpcoe)

use gugaci_global, only: jb_sys, jml, jmr, jpel, jper, jud, just, jwl, jwr, lsm_inn, norb_dz, norb_frz, norb_inn, vint_ci, voint, &
                         w0, w0_sd, w1, w1_sd
use Symmetry_Info, only: Mul
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: lra, lpcoe(norb_dz+1:norb_inn)
integer(kind=iwp) :: iwdl, iwdl1, iwdr, jmlr, kcoe, list, lmd, lmi, lmij, lmj, lr, lrd, lri, lrj, lrk, ni, nocc
real(kind=wp) :: tcoe, vlop0, vlop1, w0sd1, w0sd11, w0sd12, w0sd14, w0sd15, w0sd16, w0sd2, w0sd3, w0sd4, w0sd5, w0sd6, w0sd7, &
                 w0sd8, w0sd9, w1sd10, w1sd11, w1sd2, w1sd5, w1sd6, w1sd7, w1sd8, wl
integer(kind=iwp), external :: list3, list4

jmlr = Mul(jml,jmr)
! sd(6-1) a&r(02)-
! sd(6-2) c(22)a&r(13)-
! sd(6-4) a&r(23)c'(12)-
! sd(6-5) a&r(23)b&r(13)b^r(32)
! sd(6-8) a&r(23)b&l(32)b^l(13)
! sd(6-9) d&r&r(03)b^r(32)
! sd(6-11) d&r&l(22)b^l(13)
! sd(6-12) d&r&l(33)b^l(02)
! sd(6-13) (22)d&r&l(33)b^l(13)
! sd(6-14) d&r&l(33)c"(22)b^l(13)
! sd(6-16) d&r&l(33)b^l(23)c'(12)

! sd(6-3) a&r(13)c'(22)-
! sd(6-6) a&r(13)b&r(23)b^r(32)
! sd(6-7) a&r(13)b&l(32)b^l(23)
! sd(6-10) d&r&l(12)b^l(23)
! sd(6-15) d&r&l(33)b^l(13)c'(22)

if (jml == 1) then

  ! sd(6-1) a&r(02)-
  do lri=norb_frz+1,norb_dz
    lmi = lsm_inn(lri)
    iwdl = just(lri,lri)
    w0sd1 = w0_sd(1)
    w0sd12 = w0_sd(12)
    ni = mod(norb_dz-lri,2)
    if (ni == 1) then
      w0sd1 = -w0sd1
      w0sd12 = -w0sd12
    end if
    if (lmi == jmlr) then
      iwdr = jud(lri)
      vlop0 = w0*w0sd1
      wl = vlop0*voint(lri,lra)             !310,act_coe,610,710
      do lr=lri+1,norb_dz
        list = list3(lri,lra,lr)
        wl = wl+vlop0*(2*vint_ci(list+1)-vint_ci(list)) !  310:neoc=2,
      end do
      do lrk=norb_dz+1,lra
        list = list3(lri,lra,lrk)
        kcoe = lpcoe(lrk)
        call neoc(kcoe,nocc,tcoe)
        wl = wl+vlop0*nocc*(vint_ci(list+1)+tcoe*vint_ci(list))
      end do
      !wl = wl*vlop0
      ! sd(6-12) d&rl(33)b^l(02)
      vlop0 = w0*w0sd12
      do lrk=1,lri-1
        list = list3(lri,lra,lrk)
        wl = wl+vlop0*(vint_ci(list)-2*vint_ci(list+1))
      end do
      call prodab(3,jpel,iwdl,iwdr,jwl,jwr,wl,jper)
    end if

    ! sd(6-9) d&r&r(03)b^r(32)
    do lrd=lri+1,norb_dz
      lmd = lsm_inn(lrd)
      if (lmd /= jmr) cycle
      w0sd9 = w0_sd(9)
      ni = mod(norb_dz-lrd,2)
      if (ni == 1) w0sd9 = -w0sd9
      iwdr = jud(lrd)
      vlop0 = w0*w0sd9
      list = list3(lrd,lra,lri)
      wl = vint_ci(list)*vlop0
      call prodab(3,jpel,iwdl,iwdr,jwl,jwr,wl,jper)
    end do
  end do

end if

do lri=norb_frz+1,norb_dz
  lmi = lsm_inn(lri)
  do lrj=lri+1,norb_dz
    lmj = lsm_inn(lrj)
    lmij = Mul(lmi,lmj)
    if (lmij /= jml) cycle
    w0sd2 = w0_sd(2)
    w1sd2 = w1_sd(2)
    w0sd11 = w0_sd(11)
    w1sd11 = w1_sd(11)
    w0sd14 = w0_sd(14)
    ni = mod(norb_dz-lrj,2)
    if (ni == 1) then
      w0sd2 = -w0sd2
      w1sd2 = -w1sd2
      w0sd11 = -w0sd11
      w1sd11 = -w1sd11
      w0sd14 = -w0sd14
    end if
    iwdl = just(lri,lrj)
    iwdl1 = just(lrj,lri)
    !*******************************************************************
    ! sd(6-2) c(22)-a&r(13)-
    if (lmi == jmr) then
      iwdr = jud(lri)
      vlop0 = w0*w0sd2
      list = list3(lrj,lra,lrj)
      wl = vlop0*(voint(lrj,lra)+vint_ci(list))          !310,act_co
      do lr=lrj+1,norb_dz
        list = list3(lrj,lra,lr)
        wl = wl+vlop0*(2*vint_ci(list+1)-vint_ci(list)) !  310:neoc=
      end do
      do lrk=norb_dz+1,lra
        list = list3(lrj,lra,lrk)
        kcoe = lpcoe(lrk)
        call neoc(kcoe,nocc,tcoe)
        wl = wl+vlop0*nocc*(vint_ci(list+1)+tcoe*vint_ci(list))
      end do
      !wl = wl*vlop0
      ! sd(6-11) d&r&l(22)b^l(13)
      vlop0 = w0*w0sd11
      vlop1 = w1*w1sd11
      list = list3(lrj,lra,lri)
      wl = wl+(vlop0-vlop1)*vint_ci(list)-2*vlop0*vint_ci(list+1)
      ! sd(6-13) (22)d&r&l(33)b^l(13)
      ! sd(6-14) d&r&l(33)c"(22)b^l(13)
      vlop0 = w0*w0sd14
      do lrk=1,lrj-1
        if (lrk == lri) cycle
        list = list3(lrj,lra,lrk)
        wl = wl+vlop0*(vint_ci(list)-2*vint_ci(list+1))
      end do
      call prodab(3,jpel,iwdl,iwdr,jwl,jwr,wl,jper)
    end if
    if (jb_sys > 0) then
      if (lmj == jmr) then
        iwdr = jud(lrj)
        w0sd3 = w0_sd(3)
        w0sd15 = w0_sd(15)
        w1sd10 = w1_sd(10)
        ni = mod(norb_dz-lri,2)
        if (ni == 0) then
          w0sd3 = -w0sd3
          w0sd15 = -w0sd15
          w1sd10 = -w1sd10
        end if
        ! sd(6-3) a&r(13)-c'(22)-
        vlop0 = w0*w0sd3
        list = list3(lri,lra,lri)
        wl = voint(lri,lra)+vint_ci(list)            !310,act_coe,610
        list = list3(lri,lra,lrj)
        wl = wl+vint_ci(list+1)-vint_ci(list)  !310 c(22)neoc=1,coe=
        do lr=lri+1,norb_dz
          if (lr == lrj) cycle
          list = list3(lri,lra,lr)
          wl = wl+2*vint_ci(list+1)-vint_ci(list)     !310:neoc=2,co
        end do
        do lrk=norb_dz+1,lra
          list = list3(lri,lra,lrk)
          kcoe = lpcoe(lrk)
          call neoc(kcoe,nocc,tcoe)
          wl = wl+nocc*(vint_ci(list+1)+tcoe*vint_ci(list))
        end do
        wl = wl*vlop0
        ! sd(6-15) d&r&l(33)b^l(13)c'(22)
        vlop0 = w0*w0sd15
        do lrk=1,lri-1
          list = list3(lri,lra,lrk)
          wl = wl-vlop0*(2*vint_ci(list+1)-vint_ci(list))
        end do
        call prodab(3,jpel,iwdl1,iwdr,jwl,jwr,wl,jper)
      end if
      if ((lmij == jml) .and. (lmi == jmr)) then
        iwdr = jud(lri)
        iwdl1 = just(lrj,lri)
        w1sd10 = w1_sd(10)
        if (mod(norb_dz-lrj,2) == 1) then
          w1sd10 = -w1sd10
        end if
        ! sd(6-10) d&r&l(12)b^l(23)
        vlop1 = w1*w1sd10
        list = list3(lrj,lra,lri)
        wl = -vlop1*vint_ci(list)      !4.3
        call prodab(3,jpel,iwdl1,iwdr,jwl,jwr,wl,jper)
      end if
    end if
    ! sd(6-4) a&r(23)-c'(12)-
    if (lmj == jmr) then
      iwdr = jud(lrj)
      w0sd4 = w0_sd(4)
      w0sd16 = w0_sd(16)
      ni = mod(norb_dz-lri,2)
      if (ni == 0) w0sd4 = -w0sd4
      if (ni == 0) w0sd16 = -w0sd16
      vlop0 = w0*w0sd4
      list = list3(lri,lra,lri)
      wl = voint(lri,lra)+vint_ci(list)             !310,act_coe,610
      list = list3(lri,lra,lrj)
      wl = wl+vint_ci(list+1)-(jb_sys+2)*vint_ci(list)
      do lr=lri+1,norb_dz
        if (lr == lrj) cycle
        list = list3(lri,lra,lr)
        wl = wl+2*vint_ci(list+1)-vint_ci(list)       !310:neoc=2,co
      end do
      do lrk=norb_dz+1,lra
        list = list3(lri,lra,lrk)
        kcoe = lpcoe(lrk)
        call neoc(kcoe,nocc,tcoe)
        wl = wl+nocc*(vint_ci(list+1)+tcoe*vint_ci(list))
      end do
      wl = wl*vlop0
      ! sd(6-16) d&r&l(33)b^l(23)c'(12)
      vlop0 = w0*w0sd16
      do lrk=1,lri-1
        list = list3(lri,lra,lrk)
        wl = wl-vlop0*(2*vint_ci(list+1)-vint_ci(list))
      end do
      call prodab(3,jpel,iwdl,iwdr,jwl,jwr,wl,jper)
    end if

    ! sd(6-5) a&r(23)b&r(13)b^r(32)
    do lrd=lrj+1,norb_dz
      lmd = lsm_inn(lrd)
      if (lmd /= jmr) cycle
      iwdr = jud(lrd)
      w0sd5 = w0_sd(5)
      w1sd5 = w1_sd(5)
      ni = mod(lrj-lri+norb_dz-lrd,2)
      if (ni == 0) w0sd5 = -w0sd5
      if (ni == 0) w1sd5 = -w1sd5
      vlop0 = w0*w0sd5
      vlop1 = w1*w1sd5
      list = list4(lri,lrj,lrd,lra)
      wl = (vlop0-vlop1)*vint_ci(list+2)+(vlop0+vlop1)*vint_ci(list)     !1.3
      call prodab(3,jpel,iwdl,iwdr,jwl,jwr,wl,jper)
    end do
    if (jb_sys > 0) then
      ! sd(6-6) a&r(13)b&r(23)b^r(32)
      do lrd=lrj+1,norb_dz
        lmd = lsm_inn(lrd)
        if (lmd /= jmr) cycle
        iwdr = jud(lrd)
        w0sd6 = w0_sd(6)
        w1sd6 = w1_sd(6)
        ni = mod(lrj-lri+norb_dz-lrd,2)
        if (ni == 0) w0sd6 = -w0sd6
        if (ni == 0) w1sd6 = -w1sd6
        vlop0 = w0*w0sd6
        vlop1 = w1*w1sd6
        list = list4(lri,lrj,lrd,lra)
        wl = (vlop0-vlop1)*vint_ci(list+2)+(vlop0+vlop1)*vint_ci(list)   !1.3
        call prodab(3,jpel,iwdl1,iwdr,jwl,jwr,wl,jper)
      end do
      ! sd(6-7) a&r(13)b&l(32)b^l(23)
      do lrd=lri+1,lrj-1
        lmd = lsm_inn(lrd)
        if (lmd /= jmr) cycle
        iwdr = jud(lrd)
        w0sd7 = w0_sd(7)
        w1sd7 = w1_sd(7)
        ni = mod(lrd-lri+norb_dz-lrj,2)
        if (ni == 0) w0sd7 = -w0sd7
        if (ni == 0) w1sd7 = -w1sd7
        vlop0 = w0*w0sd7
        vlop1 = w1*w1sd7
        list = list4(lri,lrd,lrj,lra)
        wl = (vlop0-vlop1)*vint_ci(list+2)-2*vlop0*vint_ci(list+1)       !1.2
        call prodab(3,jpel,iwdl1,iwdr,jwl,jwr,wl,jper)
      end do
    end if
    ! sd(6-8) a&r(23)b&l(32)b^l(13)
    do lrd=lri+1,lrj-1
      lmd = lsm_inn(lrd)
      if (lmd /= jmr) cycle
      iwdr = jud(lrd)
      w0sd8 = w0_sd(8)
      w1sd8 = w1_sd(8)
      ni = mod(lrd-lri+norb_dz-lrj,2)
      if (ni == 0) w0sd8 = -w0sd8
      if (ni == 0) w1sd8 = -w1sd8
      vlop0 = w0*w0sd8
      vlop1 = w1*w1sd8
      list = list4(lri,lrd,lrj,lra)
      wl = (vlop0-vlop1)*vint_ci(list+2)-2*vlop0*vint_ci(list+1)         !1.2
      call prodab(3,jpel,iwdl,iwdr,jwl,jwr,wl,jper)
    end do
  end do
end do

return

end subroutine sd_head_dbl_tail_act

subroutine sdd_head_dbl_tail_act(lra,lpcoe)

use gugaci_global, only: jb_sys, jml, jmr, jpel, jper, jud, just, jwl, jwr, lsm_inn, norb_dz, norb_frz, norb_inn, vint_ci, voint, &
                         w0, w0_sd1, w1, w1_sd1
use Symmetry_Info, only: Mul
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: lra, lpcoe(norb_dz+1:norb_inn)
integer(kind=iwp) :: iwdl, iwdl1, iwdr, jmlr, kcoe, list, lmd, lmi, lmij, lmj, lr, lrd, lri, lrj, lrk, ni, nocc
real(kind=wp) :: tcoe, vlop0, vlop1, w0sd1, w0sd10, w0sd11, w0sd12, w0sd13, w0sd2, w0sd3, w0sd4, w0sd5, w0sd6, w0sd7, w0sd8, &
                 w0sd9, w1sd10, w1sd2, w1sd3, w1sd4, w1sd5, w1sd6, w1sd7, w1sd8, wl
integer(kind=iwp), external :: list3, list4

jmlr = Mul(jml,jmr)
! sd1(8-1) a&r(01)-
! sd1(8-2) c(11)a&r(23)-
! sd1(8-3) a&r(13)c'(21)-
! sd1(8-4) a&r(23)c'(11)-
! sd1(8-5) a&r(13)b&r(23)b^r(31)
! sd1(8-6) a&r(23)b&r(13)b^r(31)
! sd1(8-7) a&r(13)b&l(31)b^l(23)
! sd1(8-8) a&r(23)b&l(31)b^l(13)
! sd1(8-9) d&r&r(03)b^r(31)
! sd1(8-10) d&r&l(11)b^l(23)
! sd1(8-11) d&r&l(33)b^l(01)
! sd1(8-12) d&r&l(33)b^l(23)
! sd1(8-13) d&r&l(33)c"(13)b^l(23)
! sd1(8-14) d&r&l(33)b^l(11)c'(23)
! sd1(8-15) d&r&l(33)b^l(23)c'(11)
if (jml == 1) then
  do lri=norb_frz+1,norb_dz
    ! sd1(8-1) a&r(01)-
    lmi = lsm_inn(lri)
    iwdl = just(lri,lri)
    w0sd1 = w0_sd1(1)
    w0sd11 = w0_sd1(9)
    ni = mod(norb_dz-lri,2)
    if (ni == 1) then
      w0sd1 = -w0sd1
      w0sd11 = -w0sd11
    end if
    if (lmi == jmlr) then
      iwdr = jud(lri)
      vlop0 = w0*w0sd1
      wl = vlop0*voint(lri,lra)             !310,act_coe,610,710
      do lr=lri+1,norb_dz
        list = list3(lri,lra,lr)
        wl = wl+vlop0*(2*vint_ci(list+1)-vint_ci(list)) !  310:neoc=2,
      end do
      do lrk=norb_dz+1,lra
        list = list3(lri,lra,lrk)
        kcoe = lpcoe(lrk)
        call neoc(kcoe,nocc,tcoe)
        wl = wl+vlop0*nocc*(vint_ci(list+1)+tcoe*vint_ci(list))
      end do
      !wl = wl*vlop0
      ! sd1(8-11) d&rl(33)b^l(01)
      vlop0 = w0*w0sd11
      do lrk=1,lri-1
        list = list3(lri,lra,lrk)
        wl = wl+vlop0*(vint_ci(list)-2*vint_ci(list+1))
      end do
      call prodab(3,jpel,iwdl,iwdr,jwl,jwr,wl,jper)
    end if
    ! sd1(8-9) d&r&r(03)b^r(31)
    do lrd=lri+1,norb_dz
      lmd = lsm_inn(lrd)
      if (lmd /= jmr) cycle
      w0sd9 = w0_sd1(9)
      ni = mod(norb_dz-lrd,2)
      if (ni == 1) w0sd9 = -w0sd9
      iwdr = jud(lrd)
      vlop0 = w0*w0sd9
      list = list3(lrd,lra,lri)
      wl = vint_ci(list)*vlop0
      call prodab(3,jpel,iwdl,iwdr,jwl,jwr,wl,jper)
    end do
  end do
end if

do lri=norb_frz+1,norb_dz
  lmi = lsm_inn(lri)
  do lrj=lri+1,norb_dz
    if (lri == lrj) cycle
    lmj = lsm_inn(lrj)
    lmij = Mul(lmi,lmj)
    if (lmij /= jml) cycle
    w0sd2 = w0_sd1(2)
    w1sd2 = w1_sd1(2)
    w0sd3 = w0_sd1(3)
    w1sd3 = w1_sd1(3)
    w0sd4 = w0_sd1(4)
    w1sd4 = w1_sd1(4)
    w0sd10 = w0_sd1(10)
    w1sd10 = w1_sd1(10)
    w0sd11 = w0_sd1(11)
    w0sd12 = w0_sd1(12)
    w0sd13 = w0_sd1(13)
    ni = mod(norb_dz-lrj,2)
    if (ni == 1) then
      w0sd2 = -w0sd2
      w1sd2 = -w1sd2
      w0sd10 = -w0sd10
      w1sd10 = -w1sd10
      w0sd11 = -w0sd11
    end if
    if (mod(norb_dz-lri,2) == 1) then
      w0sd3 = -w0sd3
      w1sd3 = -w1sd3
      w0sd4 = -w0sd4
      w1sd4 = -w1sd4
      w0sd12 = -w0sd12
      w0sd13 = -w0sd13
    end if
    iwdl = just(lrj,lri)
    iwdl1 = just(lri,lrj)
    !*******************************************************************
    ! sd1(8-2) c(11)-a&r(23)-
    if (lmi == jmr) then
      iwdl = just(lrj,lri)
      iwdr = jud(lri)
      vlop0 = w0*w0sd2
      list = list3(lrj,lra,lrj)
      wl = vlop0*(voint(lrj,lra)+vint_ci(list))          !310,act_co
      do lr=lrj+1,norb_dz
        list = list3(lrj,lra,lr)
        wl = wl+vlop0*(2*vint_ci(list+1)-vint_ci(list)) !  310:neoc=
      end do
      do lrk=norb_dz+1,lra
        list = list3(lrj,lra,lrk)
        kcoe = lpcoe(lrk)
        call neoc(kcoe,nocc,tcoe)
        wl = wl+vlop0*nocc*(vint_ci(list+1)+tcoe*vint_ci(list))
      end do
      !wl = wl*vlop0
      ! sd1(8-10) d&r&l(11)b^l(23)
      vlop0 = w0*w0sd10
      vlop1 = w1*w1sd10
      list = list3(lrj,lra,lri)
      wl = wl+(vlop0-vlop1)*vint_ci(list)-2*vlop0*vint_ci(list+1)
      ! sd1(8-11) (11)d&r&l(33)b^l(23)
      ! sd1(8-11) d&r&l(33)c"(11)b^l(23)
      vlop0 = w0*w0sd11
      do lrk=1,lrj-1
        if (lrk == lri) cycle
        list = list3(lrj,lra,lrk)
        wl = wl+vlop0*(vint_ci(list)-2*vint_ci(list+1))
      end do
      call prodab(3,jpel,iwdl,iwdr,jwl,jwr,wl,jper)
    end if

    ! sd1(8-3) a&r(13)-c'(21)-
    if (lmj == jmr) then
      vlop0 = -w0*w0sd3
      list = list3(lri,lra,lri)
      wl = voint(lri,lra)+vint_ci(list)             !310,act_coe,610
      list = list3(lri,lra,lrj)
      wl = wl+vint_ci(list+1)+jb_sys*vint_ci(list)         !310
      do lr=lri+1,norb_dz
        if (lr == lrj) cycle
        list = list3(lri,lra,lr)
        wl = wl+2*vint_ci(list+1)-vint_ci(list)       !310:neoc=2,co
      end do
      do lrk=norb_dz+1,lra
        list = list3(lri,lra,lrk)
        kcoe = lpcoe(lrk)
        call neoc(kcoe,nocc,tcoe)
        wl = wl+nocc*(vint_ci(list+1)+tcoe*vint_ci(list))
      end do
      wl = wl*vlop0
      ! sd1(8-12) drl(33)-bl(13)-c'(21)-
      iwdl = just(lrj,lri)
      iwdr = jud(lrj)
      vlop0 = -w0*w0sd12
      do lrk=1,lri-1
        list = list3(lri,lra,lrk)
        wl = wl-vlop0*(2*vint_ci(list+1)-vint_ci(list))
      end do
      call prodab(3,jpel,iwdl,iwdr,jwl,jwr,wl,jper)
    end if
    ! sd1(8-4) a&r(23)-c'(11)-
    if (lmj == jmr) then
      vlop0 = -w0*w0sd4
      list = list3(lri,lra,lri)
      wl = voint(lri,lra)+vint_ci(list)           !310,act_coe,610
      list = list3(lri,lra,lrj)
      wl = wl+vint_ci(list+1)-vint_ci(list)       !310 c(11)neoc=
      do lr=lri+1,norb_dz
        if (lr == lrj) cycle
        list = list3(lri,lra,lr)
        wl = wl+2*vint_ci(list+1)-vint_ci(list)   !310:neoc=2,co
      end do
      do lrk=norb_dz+1,lra
        list = list3(lri,lra,lrk)
        kcoe = lpcoe(lrk)
        call neoc(kcoe,nocc,tcoe)
        wl = wl+nocc*(vint_ci(list+1)+tcoe*vint_ci(list))
      end do
      wl = wl*vlop0
      ! sd1(8-13) drl(33)-bl(23)-c'(11)-
      iwdl = just(lri,lrj)
      iwdr = jud(lrj)
      vlop0 = -w0*w0sd13
      do lrk=1,lri-1
        list = list3(lri,lra,lrk)
        wl = wl-vlop0*(2*vint_ci(list+1)-vint_ci(list))
      end do
      call prodab(3,jpel,iwdl,iwdr,jwl,jwr,wl,jper)
    end if
    ! sd1(8-5) a&r(13)b&r(23)b^r(31)
    do lrd=lrj+1,norb_dz
      lmd = lsm_inn(lrd)
      if (lmd /= jmr) cycle
      iwdl = just(lrj,lri)
      iwdr = jud(lrd)
      w0sd5 = w0_sd1(5)
      w1sd5 = w1_sd1(5)
      ni = mod(lrj-lri+norb_dz-lrd,2)
      if (ni == 0) w0sd5 = -w0sd5
      if (ni == 0) w1sd5 = -w1sd5
      vlop0 = w0*w0sd5
      vlop1 = w1*w1sd5
      list = list4(lri,lrj,lrd,lra)
      wl = (vlop0-vlop1)*vint_ci(list+2)+(vlop0+vlop1)*vint_ci(list)     !1.3
      call prodab(3,jpel,iwdl,iwdr,jwl,jwr,wl,jper)
      ! sd1(8-6) a&r(23)b&r(13)b^r(31)
      iwdl = just(lri,lrj)
      w0sd6 = w0_sd1(6)
      w1sd6 = w1_sd1(6)
      ni = mod(lrj-lri+norb_dz-lrd,2)
      if (ni == 0) w0sd6 = -w0sd6
      if (ni == 0) w1sd6 = -w1sd6
      vlop0 = w0*w0sd6
      vlop1 = w1*w1sd6
      list = list4(lri,lrj,lrd,lra)
      wl = (vlop0-vlop1)*vint_ci(list+2)+(vlop0+vlop1)*vint_ci(list)     !1.3
      call prodab(3,jpel,iwdl1,iwdr,jwl,jwr,wl,jper)
    end do
    ! sd1(8-7) a&r(13)b&l(31)b^l(23)
    do lrd=lri+1,lrj-1
      lmd = lsm_inn(lrd)
      if (lmd /= jmr) cycle
      iwdr = jud(lrd)
      iwdl = just(lrj,lri)
      w0sd7 = w0_sd1(7)
      w1sd7 = w1_sd1(7)
      ni = mod(lrd-lri+norb_dz-lrj,2)
      if (ni == 0) w0sd7 = -w0sd7
      if (ni == 0) w1sd7 = -w1sd7
      vlop0 = w0*w0sd7
      vlop1 = w1*w1sd7
      list = list4(lri,lrd,lrj,lra)
      wl = (vlop0-vlop1)*vint_ci(list+2)-2*vlop0*vint_ci(list+1)         !1.2
      call prodab(3,jpel,iwdl,iwdr,jwl,jwr,wl,jper)
      ! sd1(8-8) a&r(23)b&l(31)b^l(13)
      iwdl = just(lri,lrj)
      w0sd8 = w0_sd1(8)
      w1sd8 = w1_sd1(8)
      ni = mod(lrd-lri+norb_dz-lrj,2)
      if (ni == 0) w0sd8 = -w0sd8
      if (ni == 0) w1sd8 = -w1sd8
      vlop0 = w0*w0sd8
      vlop1 = w1*w1sd8
      list = list4(lri,lrd,lrj,lra)
      wl = (vlop0-vlop1)*vint_ci(list+2)-2*vlop0*vint_ci(list+1)         !1.2
      call prodab(3,jpel,iwdl1,iwdr,jwl,jwr,wl,jper)
    end do
  end do
end do

return

end subroutine sdd_head_dbl_tail_act

!subroutine td_head_dbl_tail_act(lra,lpcoe)
!
!use gugaci_global, only: jml, jmr, jpel, jper, jud, just, jwl, jwr, lsm_inn, norb_dz, norb_frz, norb_inn, vint_ci, voint, w0, &
!                         w0_td, w1, w1_td
!use Symmetry_Info, only: Mul
!use Definitions, only: wp, iwp
!
!implicit none
!integer(kind=iwp), intent(in) :: lra, lpcoe(norb_dz+1:norb_inn)
!integer(kind=iwp) :: iwdl, iwdr, jmlr, kcoe, list, lmd, lmi, lmij, lmj, lr, lrd, lri, lrj, lrk, ni, nocc
!real(kind=wp) :: tcoe, vlop0, vlop1, w0td1, w0td2, w0td3, w0td4, w0td5, w1td2, w1td3, wl
!integer(kind=iwp), external :: list3, list4
!
!jmlr = Mul(jml,jmr)
!! td(13-1) (22)a&(23)
!! td(13-1) a&(23)c'(22)
!! td(13-5) (22)d&&l(33)b^l(23)
!do lri=norb_frz+1,norb_dz
!  lmi = lsm_inn(lri)
!  if (lmi /= jmlr) cycle
!  w0td1 = w0_td(1)
!  w0td4 = w0_td(4)
!  w0td5 = w0_td(5)
!  ni = mod(norb_dz-lri,2)
!  if (ni == 1) w0td1 = -w0td1
!  if (ni == 1) w0td4 = -w0td4
!  if (ni == 1) w0td5 = -w0td5
!
!  ! td(13-1) a&(23)c'(22)
!  do lrd=lri+1,norb_dz
!    lmd = lsm_inn(lrd)
!    if (lmd /= jmr) cycle
!    iwdl = just(lri,lrd)
!    iwdr = jud(lrd)
!    vlop0 = -w0*w0td1
!    list = list3(lri,lra,lri)
!    wl = voint(lri,lra)+vint_ci(list)        !310,act_coe,610,7
!    list = list3(lri,lra,lrd)
!    wl = wl+vint_ci(list+1)                  !310 c(22) coe
!    do lr=lri+1,norb_dz
!      if (lr == lrd) cycle
!      list = list3(lri,lra,lr)
!      wl = wl+2*vint_ci(list+1)-vint_ci(list)       !310:neoc=2,coe=
!    end do
!    do lrk=norb_dz+1,lra
!      list = list3(lri,lra,lrk)
!      kcoe = lpcoe(lrk)
!      call neoc(kcoe,nocc,tcoe)
!      wl = wl+nocc*(vint_ci(list+1)+tcoe*vint_ci(list))
!    end do
!    wl = wl*vlop0
!    ! td(13-5) d&rl(33)b^l(23)c'(22)          !cc (22)d&rl(33)b^l(23)???
!    vlop0 = -w0*w0td5
!    do lrk=1,lri-1
!      list = list3(lri,lra,lrk)
!      wl = wl-vlop0*(2*vint_ci(list+1)-vint_ci(list))
!    end do
!    call prodab(3,jpel,iwdl,iwdr,jwl,jwr,wl,jper)
!  end do
!  !---------------------------------------------------------------------
!  do lrd=norb_frz+1,lri-1
!    lmd = lsm_inn(lrd)
!    if (lmd /= jmr) cycle
!    iwdl = just(lrd,lri)
!    iwdr = jud(lrd)
!    ! td(13-1) (22)a&(23)
!    vlop0 = w0*w0td1
!    list = list3(lri,lra,lri)
!    wl = vlop0*(voint(lri,lra)+vint_ci(list))             !310,act_c
!    do lr=lri+1,norb_dz
!      list = list3(lri,lra,lr)
!      wl = wl+vlop0*(2*vint_ci(list+1)-vint_ci(list))     !310:neoc=2,
!    end do
!    do lrk=norb_dz+1,lra
!      list = list3(lri,lra,lrk)
!      kcoe = lpcoe(lrk)
!      call neoc(kcoe,nocc,tcoe)
!      wl = wl+vlop0*nocc*(vint_ci(list+1)+tcoe*vint_ci(list))
!    end do
!    !wl = wl*vlop0
!    ! td(13-4) d&r&l(22)b^l(23)
!    vlop0 = w0*w0td4
!    vlop1 = w1*w0td4
!    list = list3(lri,lra,lrd)
!    wl = wl+(vlop0-vlop1)*vint_ci(list)-2*vlop0*vint_ci(list+1)
!    ! td(13-5) d&rl(33)c"(22)b^l(23)
!    vlop0 = w0*w0td5
!    do lrk=1,lri-1
!      if (lrk == lrd) cycle
!      list = list3(lri,lra,lrk)
!      wl = wl+vlop0*(vint_ci(list)-2*vint_ci(list+1))      !4.3
!    end do
!    call prodab(3,jpel,iwdl,iwdr,jwl,jwr,wl,jper)
!  end do
!end do
!do lri=norb_frz+1,norb_dz-1
!  lmi = lsm_inn(lri)
!  do lrj=lri+1,norb_dz
!    lmj = lsm_inn(lrj)
!    lmij = Mul(lmi,lmj)
!    if (lmij /= jml) cycle
!    iwdl = just(lri,lrj)
!
!    ! td(13-2) a&(23)b&r(23)b^r(32)
!    do lrd=lrj+1,norb_dz
!      lmd = lsm_inn(lrd)
!      if (lmd /= jmr) cycle
!      w0td2 = w0_td(2)
!      w1td2 = w1_td(2)
!      ni = mod(lrj-lri+norb_dz-lrd,2)
!      if (ni == 0) w0td2 = -w0td2
!      if (ni == 0) w1td2 = -w1td2
!
!      iwdr = jud(lrd)
!      vlop0 = w0*w0td2
!      vlop1 = w1*w1td2
!      list = list4(lri,lrj,lrd,lra)
!      wl = vlop0*(vint_ci(list+2)+vint_ci(list))-vlop1*(vint_ci(list+2)-vint_ci(list)) !1.3
!      call prodab(3,jpel,iwdl,iwdr,jwl,jwr,wl,jper)
!    end do
!    ! td(13-3) a&(23)b&l(32)b^l(23)
!    do lrd=lri+1,lrj-1
!      lmd = lsm_inn(lrd)
!      if (lmd /= jmr) cycle
!      iwdr = jud(lrd)
!      w0td3 = w0_td(3)
!      w1td3 = w1_td(3)
!      ni = mod(lrd-lri+norb_dz-lrj,2)
!      if (ni == 0) w0td3 = -w0td3
!      if (ni == 0) w1td3 = -w1td3
!      vlop0 = w0*w0td3                !d6-8
!      vlop1 = w1*w1td3
!      list = list4(lri,lrd,lrj,lra)
!      wl = vlop0*(vint_ci(list+2)-2*vint_ci(list+1))-vlop1*vint_ci(list+2)    !1.2
!      call prodab(3,jpel,iwdl,iwdr,jwl,jwr,wl,jper)
!    end do
!  end do
!end do
!
!return
!
!end subroutine td_head_dbl_tail_act

subroutine ttdd_head_dbl_tail_act(lra,lpcoe)

use gugaci_global, only: jml, jmr, jpel, jper, jud, just, jwl, jwr, lsm_inn, norb_dz, norb_frz, norb_inn, vint_ci, voint, w0, &
                         w0_t1d1, w1, w1_t1d1
use Symmetry_Info, only: Mul
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: lra, lpcoe(norb_dz+1:norb_inn)
integer(kind=iwp) :: iwdl, iwdr, jmlr, kcoe, list, lmd, lmi, lmij, lmj, lr, lrd, lri, lrj, lrk, ni, nocc
real(kind=wp) :: tcoe, vlop0, vlop1, w0td1, w0td2, w0td3, w0td4, w0td5, w1td2, w1td3, wl
integer(kind=iwp), external :: list3, list4

jmlr = Mul(jml,jmr)
! t1d1(15-1) (11)a&(13)
! t1d1(15-1) a&(13)c'(11)
! t1d1(15-5) (11)d&&l(33)b^l(13)
do lri=norb_frz+1,norb_dz
  lmi = lsm_inn(lri)
  if (lmi /= jmlr) cycle
  w0td1 = w0_t1d1(1)
  w0td4 = w0_t1d1(4)
  w0td5 = w0_t1d1(5)
  ni = mod(norb_dz-lri,2)
  if (ni == 1) then
    w0td1 = -w0td1
    w0td4 = -w0td4
    w0td5 = -w0td5
  end if
  ! t1d1(15-2) a&(13)c'(11)
  do lrd=lri+1,norb_dz
    lmd = lsm_inn(lrd)
    if (lmd /= jmr) cycle
    iwdl = just(lri,lrd)
    iwdr = jud(lrd)
    vlop0 = -w0*w0td1
    list = list3(lri,lra,lri)
    wl = voint(lri,lra)+vint_ci(list)           !310,act_coe,610,7
    list = list3(lri,lra,lrd)
    wl = wl+vint_ci(list+1)
    do lr=lri+1,norb_dz
      if (lr == lrd) cycle
      list = list3(lri,lra,lr)
      wl = wl+2*vint_ci(list+1)-vint_ci(list)       !310:neoc=2,coe=
    end do
    do lrk=norb_dz+1,lra
      list = list3(lri,lra,lrk)
      kcoe = lpcoe(lrk)
      call neoc(kcoe,nocc,tcoe)
      wl = wl+nocc*(vint_ci(list+1)+tcoe*vint_ci(list))
    end do
    wl = wl*vlop0
    ! t1d1(15-5) d&rl(33)b^l(13)c'(11)
    vlop0 = -w0*w0td5
    do lrk=1,lri-1
      list = list3(lri,lra,lrk)
      wl = wl-vlop0*(2*vint_ci(list+1)-vint_ci(list))
    end do
    call prodab(3,jpel,iwdl,iwdr,jwl,jwr,wl,jper)
  end do
  !---------------------------------------------------------------------
  do lrd=norb_frz+1,lri-1
    lmd = lsm_inn(lrd)
    if (lmd /= jmr) cycle
    iwdl = just(lrd,lri)
    iwdr = jud(lrd)
    ! t1d1(15-1) (11)a&(13)
    vlop0 = w0*w0td1
    list = list3(lri,lra,lri)
    wl = vlop0*(voint(lri,lra)+vint_ci(list))             !310,act_c
    do lr=lri+1,norb_dz
      list = list3(lri,lra,lr)
      wl = wl+vlop0*(2*vint_ci(list+1)-vint_ci(list))     !310:neoc=2,
    end do
    do lrk=norb_dz+1,lra
      list = list3(lri,lra,lrk)
      kcoe = lpcoe(lrk)
      call neoc(kcoe,nocc,tcoe)
      wl = wl+vlop0*nocc*(vint_ci(list+1)+tcoe*vint_ci(list))
    end do
    !wl = wl*vlop0
    ! t1d1(15-4) d&r&l(11)b^l(13)
    vlop0 = w0*w0td4
    vlop1 = w1*w0td4
    list = list3(lri,lra,lrd)
    wl = wl+(vlop0-vlop1)*vint_ci(list)-2*vlop0*vint_ci(list+1)
    ! t1d1(15-5) d&rl(33)c"(11)b^l(13)
    vlop0 = w0*w0td5
    do lrk=1,lri-1
      if (lrk == lrd) cycle
      list = list3(lri,lra,lrk)
      wl = wl+vlop0*(vint_ci(list)-2*vint_ci(list+1))      !4.3
    end do
    call prodab(3,jpel,iwdl,iwdr,jwl,jwr,wl,jper)
  end do
end do
do lri=norb_frz+1,norb_dz-1
  lmi = lsm_inn(lri)
  do lrj=lri+1,norb_dz
    lmj = lsm_inn(lrj)
    lmij = Mul(lmi,lmj)
    if (lmij /= jml) cycle
    iwdl = just(lri,lrj)

    ! t1d1(15-2) a&(13)b&r(13)b^r(31)
    do lrd=lrj+1,norb_dz
      lmd = lsm_inn(lrd)
      if (lmd /= jmr) cycle
      w0td2 = w0_t1d1(2)
      w1td2 = w1_t1d1(2)
      ni = mod(lrj-lri+norb_dz-lrd,2)
      if (ni == 0) w0td2 = -w0td2
      if (ni == 0) w1td2 = -w1td2

      iwdr = jud(lrd)
      vlop0 = w0*w0td2
      vlop1 = w1*w1td2
      list = list4(lri,lrj,lrd,lra)
      wl = vlop0*(vint_ci(list+2)+vint_ci(list))-vlop1*(vint_ci(list+2)-vint_ci(list))   !1.3
      call prodab(3,jpel,iwdl,iwdr,jwl,jwr,wl,jper)
    end do
    ! t1d1(15-3) a&(13)b&l(31)b^l(13)
    do lrd=lri+1,lrj-1
      lmd = lsm_inn(lrd)
      if (lmd /= jmr) cycle
      iwdr = jud(lrd)
      w0td3 = w0_t1d1(3)
      w1td3 = w1_t1d1(3)
      ni = mod(lrd-lri+norb_dz-lrj,2)
      if (ni == 0) w0td3 = -w0td3
      if (ni == 0) w1td3 = -w1td3
      vlop0 = w0*w0td3                !d6-8
      vlop1 = w1*w1td3
      list = list4(lri,lrd,lrj,lra)
      wl = vlop0*(vint_ci(list+2)-2*vint_ci(list+1))-vlop1*vint_ci(list+2)   !1.2
      call prodab(3,jpel,iwdl,iwdr,jwl,jwr,wl,jper)
    end do
  end do
end do

return

end subroutine ttdd_head_dbl_tail_act

subroutine sv_arbr_act_c_ext_stv_sgt0(lin)

use gugaci_global, only: ipae, ipael, jml, jmr, jpad, jpadl, just, lp_lwei, lp_rwei, lpnew_lwei, lpnew_rwei, lsm_inn, mhlp, mtype, &
                         norb_dz, norb_frz, vplp_w0, vplp_w1, vplpnew_w0, vplpnew_w1, w0_sv, w1_sv
use Symmetry_Info, only: Mul
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: lin
integer(kind=iwp) :: iwal, iwar, iwdl, iwdr, jmlr, lmi, lmij, lmj, lri, lrj, mpl, ni
real(kind=wp) :: w0sv1, w1sv1
integer(kind=iwp), external :: iwalk_ad

! sv(10-1) ar(13)br(23)  act -c"-  tv_ext -br-ar
jmlr = Mul(jml,jmr)
do lri=norb_frz+1,norb_dz
  lmi = lsm_inn(lri)
  do lrj=lri+1,norb_dz
    lmj = lsm_inn(lrj)
    lmij = Mul(lmi,lmj)
    if (lmij /= jmlr) cycle
    !-------------------------------------------------------------------
    w0sv1 = w0_sv(1)
    w1sv1 = w1_sv(1)
    ni = mod(lrj-lri,2)
    if (ni == 0) then
      w0sv1 = -w0sv1
      w1sv1 = -w1sv1
    end if
    iwdl = just(lrj,lri)
    iwdr = 0
    do mpl=1,mhlp
      iwal = lpnew_lwei(mpl)
      iwar = lpnew_rwei(mpl)
      lp_lwei(mpl) = iwalk_ad(jpadl,ipael,iwal,iwdl)
      lp_rwei(mpl) = iwalk_ad(jpad,ipae,iwar,iwdr)
    end do
    do mpl=1,mtype
      vplp_w0(mpl) = vplpnew_w0(mpl)*w0sv1
      vplp_w1(mpl) = vplpnew_w1(mpl)*w1sv1
    end do
    if (lin == 10) then
      call ar_br_sv_ext_br_ar(lri,lrj)
    end if
    if (lin == 17) then
      call ar_br_tv_ext_br_ar(lri,lrj)
    end if
  end do
end do

return

end subroutine sv_arbr_act_c_ext_stv_sgt0
