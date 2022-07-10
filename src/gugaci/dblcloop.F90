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

! complete double occupied space loops

subroutine dbl_space_loop()

use gugaci_global, only: norb_dbl

implicit none

if (norb_dbl == 0) return
call dbl_space_loop_ijkk_sgezero()
call dbl_space_loop_ijkl_sgezero()

return

end subroutine dbl_space_loop

subroutine dbl_space_loop_ijkk_sgezero()

use gugaci_global, only: jb_sys, jud, just, lsm_inn, norb_dz, norb_frz, ns_sm, vint_ci, voint
use Symmetry_Info, only: Mul
use Constants, only: Zero, One, Half
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp) :: im, imi, imm, iwdl, iwdr, iwld, iwls, iwlt, iwrd, iwrs, iwrt, jpat, jpdd, jpds, jpds0, jpdt, jpdt1, kij, &
                     list, lmi, lmij, lmj, lr, lri, lrj, lrm, mij, ni
real(kind=wp) :: db, vl0_2, vl_0, vl_1, vlp_1, vls0, vls0_2, vls1, vls10_2, vls10_2b, vls1_a, vls_c, vlt0, vlt1, wl, wl0, wl10, &
                 wls, wls0, wls0_1, wls0_2, wls1, wls2, wls_a, wls_b, wlt, wlt0, wltmp
integer(kind=iwp), external :: list3

!==============================  g1,2,4,6,7,8 ==========================
!zz = ' doub_800_v'
wls1 = Zero
db = jb_sys
im = ns_sm
jpds0 = im+17
do lri=norb_frz+1,norb_dz-1
  !mi = lsm_inn(lri)
  iwdl = just(lri,lri)         ! drr(03)-drr(30)
  do lrj=lri+1,norb_dz
    wl = voint(lrj,lri)        ! vdint(itailorb,iheadorb)
    iwdr = just(lrj,lrj)
    call prodab(1,0,jpds0,iwdl,iwdr,0,wl,0)
  end do
end do
! bbs debug 20090722 jb_sys >0
if (jb_sys > 0) then
  ! drl(12)-drl(21)
  vlp_1 = -sqrt(db*(db+2))/(db+1)
  do lri=norb_frz+1,norb_dz
    lmi = lsm_inn(lri)
    do lrj=lri+1,norb_dz
      lmj = lsm_inn(lrj)
      lmij = Mul(lmi,lmj)
      jpds = 17+Mul(lmij,ns_sm)
      iwdl = just(lrj,lri)
      iwdr = just(lri,lrj)
      wl = -vlp_1*voint(lrj,lri)
      !write(u6,*) 'lri,lrj',lri,lrj,jpds,iwdl,iwdr
      call prodab(1,0,jpds,iwdl,iwdr,0,wl,0)
    end do
  end do
end if

do lri=norb_frz+1,norb_dz-1           !frz
  imi = lsm_inn(lri)
  !n2 = ngw2(lri-2)
  do lrj=lri+1,norb_dz
    mij = Mul(imi,lsm_inn(lrj))
    if (mij /= 1) cycle
    ni = mod(lrj-lri,2)
    !=========== down comm for 2 4 =====================================
    vl_0 = sqrt((db+2)/(db+1))
    vl_1 = sqrt(db/(db+1))
    if (ni == 0) vl_0 = -vl_0
    if (ni == 0) vl_1 = -vl_1
    vl0_2 = One
    vls0_2 = 1/(db+1)
    vls10_2 = sqrt(db*(db+2))/(db+1)
    !if (ni == 1) vl0_2 = -vl0_2
    vls0 = -Half
    vls1 = (db+3)/(2*db+2)
    vls1_a = (db-1)/(2*db+2)
    vls_c = db/(db+1)
    vlt0 = -Half
    vlt1 = -Half
    vls10_2b = sqrt(db*(db+2))/(db+1)
    if (ni == 1) then
      vl0_2 = -vl0_2
      vls0_2 = -vls0_2
      vls10_2 = -vls10_2
      vls_c = -vls_c
    end if
    if (ni == 0) then
      vls0 = -vls0
      vls1 = -vls1
      vls1_a = -vls1_a
      vls10_2b = -vls10_2b
      vlt0 = -vlt0
      vlt1 = -vlt1
    end if
    !wl0 = Zero
    !wl10 = Zero
    ! for 2,4 the common
    ! ar(23)-c"(33)-ar(10) ar(13)-c"(33)-ar(20)
    ! ar(23)-drl(33)-ar(10)  ar(13)-drl(33)-ar(20)    310
    wltmp = Zero
    do lr=lri+1,lrj-1
      list = list3(lri,lrj,lr)
      wltmp = wltmp+2*vint_ci(list+1)-vint_ci(list)
      !wl0 = wl0+vlt0*(2*vint_ci(list+1)-vint_ci(list))
      !wl10 = wl10+vlt1*(2*vint_ci(list+1)-vint_ci(list))
      !wl0 = wl0+vl0_1*(2*vint_ci(list+1)-vint_ci(list)) !310:neoc=2,
    end do
    wl0 = vl_0*wltmp
    wl10 = vl_1*wltmp
    ! neoc(k)*wg43*(vint(list+2)+coe(k)*vint(list+1))
    ! drl(33)-bl(23)-ar(10)
    ! drl(33)-bl(13)-ar(20)
    wltmp = Zero
    do lr=1,lri-1
      list = list3(lri,lrj,lr)
      wltmp = wltmp+2*vint_ci(list+1)-vint_ci(list) ! 430 w0=-vl0 w1
    end do
    wl0 = wl0+vl_0*wltmp
    wl10 = wl10+vl_1*wltmp
    ! ar(23)-bl(10)-drl(33)
    ! ar(13)-bl(20)-drl(33)
    wltmp = Zero
    do lr=lrj+1,norb_dz
      list = list3(lri,lrj,lr)
      wltmp = wltmp+2*vint_ci(list+1)-vint_ci(list)  ! 220 w0=-vl0 w
    end do
    wl0 = wl0+vl_0*wltmp
    wl10 = wl10+vl_1*wltmp
    !=========== start comm for 2 4 ====================================
    do lrm=norb_frz+1,norb_dz        !ic=1,norb_act   !frz
      imm = lsm_inn(lrm)
      im = Mul(imm,imi)
      im = Mul(im,ns_sm)
      kij = 0
      if (lrm == lrj) kij = 2
      if (lrm == lri) kij = 4
      if (lrm < lri) kij = 7
      if (lrm > lrj) kij = 6
      if ((lrm > lri) .and. (lrm < lrj)) kij = 8
      select case (kij)
        case (2)
          ! ss(1-c1)  ar(23)-ar(10)-    ar(13)-ar(20)
          iwdl = just(lri,lrj)
          iwdr = just(lrj,lrj)
          list = list3(lri,lrj,lri)
          wl = wl0+vl_0*(voint(lri,lrj)+vint_ci(list))
          call prodab(1,0,jpds0,iwdl,iwdr,0,wl,0)
          if (jb_sys > 0) then
            wl = wl10+vl_1*(voint(lri,lrj)+vint_ci(list))
            iwdl = just(lrj,lri)
            call prodab(1,0,jpds0,iwdl,iwdr,0,wl,0)
          end if

        case (4)
          ! ss(1-c2)  ar(02)-ar(31)    ar(01)-ar(32)
          iwdl = just(lri,lri)
          iwdr = just(lri,lrj)
          list = list3(lri,lrj,lrj)
          wl = wl0+vl_0*(voint(lri,lrj)+vint_ci(list))
          call prodab(1,0,jpds0,iwdl,iwdr,0,wl,0)
          if (jb_sys > 0) then
            iwdr = just(lrj,lri)
            wl = wl10+vl_1*(voint(lri,lrj)+vint_ci(list))
            call prodab(1,0,jpds0,iwdl,iwdr,0,wl,0)
          end if

        case (6)
          ! ar(23)-br(13)-drr(30)  w0=-sqrt((db+2)/(db+1)) w1=0
          ! ar(13)-br(23)-drr(30)  w0=-sqrt(db/(db+1)) w1=0
          iwdl = just(lri,lrj)
          iwdr = just(lrm,lrm)
          list = list3(lri,lrj,lrm)
          wl = -vl_0*vint_ci(list)
          call prodab(1,0,jpds0,iwdl,iwdr,0,wl,0)
          if (jb_sys > 0) then
            iwdl = just(lrj,lri)
            wl = -vl_1*vint_ci(list)
            call prodab(1,0,jpds0,iwdl,iwdr,0,wl,0)
          end if

        case (7)
          ! drr(03)br(32)br(31)   w0=-sqrt((db+2)/(db+1)) w1=0
          ! drr(03)br(31)br(32)   w0=-sqrt(db/(db+1)) w1=0
          iwdl = just(lrm,lrm)
          iwdr = just(lri,lrj)
          list = list3(lri,lrj,lrm)
          wl = -vl_0*vint_ci(list)
          call prodab(1,0,jpds0,iwdl,iwdr,0,wl,0)
          if (jb_sys > 0) then
            iwdr = just(lrj,lri)
            iwdl = just(lrm,lrm)
            wl = -vl_1*vint_ci(list)
            call prodab(1,0,jpds0,iwdl,iwdr,0,wl,0)
          end if

        case (8)
          ! ar(23)-drl(30)-al(13)   w0=sqrt((db+2)/(db+1)) w1=0
          ! ar(13)-drl(30)-al(23)   w0=sqrt(db/(db+1)) w1=0
          iwdl = just(lri,lrj)
          iwdr = just(lrm,lrm)
          list = list3(lri,lrj,lrm)
          wl = -vl_0*vint_ci(list)
          call prodab(1,0,jpds0,iwdl,iwdr,0,wl,0)
          if (jb_sys > 0) then
            ! ar(13)-drl(30)-al(23)   w0=sqrt(db/(db+1)) w1=0
            iwdl = just(lrj,lri)
            wl = -vl_1*vint_ci(list)
            call prodab(1,0,jpds0,iwdl,iwdr,0,wl,0)
          end if
          !========= start  d9(ss) d35(tt) =============================
          jpds = 17+im
          jpdt = 9+im
          jpdt1 = jpdt+24
          iwls = just(lri,lrm)
          iwrs = just(lrm,lrj)
          iwlt = iwls
          iwrt = iwrs
          ! t  ar(23)-c'(22)-cw(33)-ar(32)   w0=1 w1=0
          ! t  ar(13)-c'(11)-cw(33)-ar(31)   w0=1 w1=0
          ! s  ar(23)-c'(12)-cw(33)-ar(31)   w0=-1/(db+1) w1=0
          ! s  ar(23)-c'(11)-cw(33)-ar(32)   w0=-sqrt(db*(db+2))/(db+1) w1=0
          ! s  ar(13)-c'(21)-cw(33)-ar(32)   w0=1/(db+1) w1=0
          ! s  ar(13)-c'(22)-cw(33)-ar(31)   w0=-sqrt(db*(db+2))/(db+1) w1=0
          wlt0 = Zero
          wls0 = Zero
          wls0_1 = Zero                   ! bbs_tmp
          wls0_2 = Zero
          wltmp = Zero
          do lr=lrm+1,lrj-1
            list = list3(lri,lrj,lr)
            wltmp = wltmp+2*vint_ci(list+1)-vint_ci(list)
            !wlt0 = wlt0+vl0_2*(2*vint_ci(list+1)-vint_ci(list))
          end do
          wlt0 = wlt0+vl0_2*wltmp
          wls0 = wls0-vls0_2*wltmp
          wls0_1 = wls0_1-vls10_2*wltmp
          wls0_2 = wls0_2+vls0_2*wltmp
          !wls2 = wls0
          ! t  ar(23)-cw(33)-c'(22)-ar(32)      w0=1 w1=0
          ! t  ar(13)-cw(33)-c'(11)-ar(31)      w0=1 w1=0
          ! s  ar(23)-cw(33)-c'(12)-ar(31)      w0=-1/(db+1) w1=0
          ! s  ar(23)-cw(33)-c'(11)-ar(32)      w0=-sqrt(db*(db+2))/(db+1) w1=0
          ! s  ar(13)-cw(33)-c'(21)-ar(32)      w0=1/(db+1) w1=0
          ! s  ar(13)-cw(33)-c'(22)-ar(31)      w0=-sqrt(db*(db+2))/(db+1) w1=0
          wltmp = Zero
          do lr=lri+1,lrm-1
            list = list3(lri,lrj,lr)
            wltmp = wltmp+2*vint_ci(list+1)-vint_ci(list)
            !wlt0 = wlt0+vl0_2*(2*vint_ci(list+1)-vint_ci(list))
          end do
          wlt0 = wlt0+vl0_2*wltmp
          wls0 = wls0-vls0_2*wltmp
          wls0_1 = wls0_1-vls10_2*wltmp
          wls0_2 = wls0_2+vls0_2*wltmp
          ! t  drl(33)-bl(23)-c'(22)-ar(32) w0=-1 w1=0
          ! t  drl(33)-bl(13)-c'(11)-ar(31) w0=-1 w1=0
          ! s  drl(33)-bl(23)-c'(12)-ar(31) w0=1/(db+1)  w1=0
          ! s  drl(33)-bl(23)-c'(11)-ar(32) w0=sqrt(db*(db+2))/(db+1) w1=0
          ! s  drl(33)-bl(13)-c'(21)-ar(32) w0=-1/(db+1)  w1=0
          ! s  drl(33)-bl(13)-c'(22)-ar(31) w0=sqrt(db*(db+2))/(db+1)   w1=0
          wltmp = Zero
          do lr=1,lri-1
            list = list3(lri,lrj,lr)
            wltmp = wltmp+(vint_ci(list)-2*vint_ci(list+1))
          end do
          wlt0 = wlt0-vl0_2*wltmp                  !bbs_tmp
          wls0 = wls0+vls0_2*wltmp
          wls0_1 = wls0_1+vls10_2*wltmp
          wls0_2 = wls0_2-vls0_2*wltmp
          ! t  ar(23)-c'(22)-bl(32)-drl(33) w0=-1 w1=0
          ! t  ar(13)-c'(11)-bl(31)-drl(33) w0=-1 w1=0
          ! s  ar(23)-c'(12)-br(31)-drl(33) w0=1/(db+1)  w1=0
          ! s  ar(23)-c'(11)-br(32)-drl(33) w0=sqrt(db*(db+2))/(db+1) w1=0
          ! s  ar(13)-c'(21)-br(32)-drl(33) w0=-1/(db+1)  w1=0
          ! s  ar(13)-c'(22)-br(31)-drl(33) w0=sqrt(db*(db+2))/(db+1) w1=0
          wltmp = Zero
          do lr=lrj+1,norb_dz
            list = list3(lri,lrj,lr)
            wltmp = wltmp+(vint_ci(list)-2*vint_ci(list+1))
          end do
          !wlt0 = wlt0+vl0_2*wltmp
          !wls0 = wls0-vls0_2*wltmp
          !wls0_1 = wls0_1+vls10_2*wltmp
          !wls0_2 = wls0_2+vls0_2*wltmp
          wlt0 = wlt0-vl0_2*wltmp
          wls0 = wls0+vls0_2*wltmp
          wls0_1 = wls0_1+vls10_2*wltmp
          wls0_2 = wls0_2-vls0_2*wltmp
          !wls0 = -wlt0
          !========= up comm for 9 35 ==================================
          ! t  ar(23)-c'(22)-arw(32) w0=1 w1=0
          ! t  ar(13)-c'(11)-arw(31) w0=1 w1=0
          ! s  ar(23)-c'(12)-arw(31) w0=-1/(db+1)  w1=0
          ! s  ar(23)-c'(11)-arw(31) w0=-sqrt(db*(db+2))/(db+1) w1=0
          ! s  ar(13)-c'(21)-arw(32) w0=1/(db+1)  w1=0
          ! s  ar(13)-c'(22)-arw(31) w0=-sqrt(db*(db+2))/(db+1) w1=0
          !list = list3(lri,lrj,lrj)
          !wlt = wlt0+vl0_2*vint_ci(list)
          !wls = wls0-vl0_2*vint_ci(list)
          list = list3(lri,lrj,lrj)
          wlt = wlt0+vl0_2*vint_ci(list)
          wls = wls0-vls0_2*vint_ci(list)
          wls1 = wls0_1-vls10_2*vint_ci(list)
          wls2 = wls0_2+vls0_2*vint_ci(list)
          !wl_bbs = -vls0_2*vint_ci(list)
          ! t  ar(23)-cw(22)-ar(32)
          ! t  ar(13)-cw(11)-ar(31)
          ! s  ar(23)-cw(12)-ar(31)    w0=-1/(db+1)    ar(12)-drr(22)-ar(32)
          ! s  ar(23)-cw(11)-ar(32)    w0=-sqrt(db*(db+1))/(db+1)
          ! s  ar(13)-cw(21)-ar(32)    w0=1/(db+1)
          ! s  ar(13)-cw(22)-ar(31)    w0=-sqrt(db*(db+1))/(db+1)
          list = list3(lri,lrj,lrm)
          wlt = wlt+vl0_2*vint_ci(list+1)
          wls = wls-vls0_2*(vint_ci(list+1)-(db+2)*vint_ci(list)) !coe=
          wls1 = wls1-vls10_2*(vint_ci(list+1)-vint_ci(list))     !coe=-1
          wls2 = wls2+vls0_2*(vint_ci(list+1)+db*vint_ci(list))   !coe=d
          !wls2 = wls2+vls0_2*(vint_ci(list+1)-2*vint_ci(list))
          !wl_bbs = wl_bbs-vls0_2*(vint_ci(list+1)-2*vint_ci(list))
          ! t  arw(23)-c'(22)-ar(32)
          ! t  arw(13)-c'(11)-ar(31)
          ! s  arw(23)-c'(12)-ar(31)
          ! s  arw(23)-c'(11)-ar(32)
          ! s  arw(13)-c'(21)-ar(32)
          ! s  arw(13)-c'(22)-ar(31)
          ! t  ar(23)-c'(22)-ar(32)   w0=1
          ! t  ar(13)-c'(11)-ar(31)   w0=1
          ! s  ar(23)-c'(12)-ar(31)   w0=-1/(db+1)  w1=0
          ! s  ar(23)-c'(11)-ar(32)   w0=-sqrt(db*(db+2))/(db+1) w1=0
          ! s  ar(13)-c'(21)-ar(32)   w0=1/(db+1)  w1=0
          ! s  ar(13)-c'(22)-ar(31)   w0=-sqrt(db*(db+2))/(db+1) w1=0
          list = list3(lri,lrj,lri)
          wlt = wlt+vl0_2*(voint(lri,lrj)+vint_ci(list))
          wls = wls-vls0_2*(voint(lri,lrj)+vint_ci(list))
          wls1 = wls1-vls10_2*(voint(lri,lrj)+vint_ci(list))
          wls2 = wls2+vls0_2*(voint(lri,lrj)+vint_ci(list))
          !wl_bbs = wl_bbs-vls0_2*(voint(lri,lrj)+vint_ci(list))

          ! ar(23)-drr(12)-ar(31)  w0=(db+2)/(db+1)   ???? complete
          call prodab(1,0,jpds,iwls,iwrs,0,wls,0)
          call prodab(1,0,jpdt,iwlt,iwrt,0,wlt,0)
          if (jb_sys > 0) then
            iwls = just(lri,lrm)
            iwrs = just(lrj,lrm)
            call prodab(1,0,jpds,iwls,iwrs,0,wls1,0)
            iwls = just(lrm,lri)
            iwrs = just(lrj,lrm)
            call prodab(1,0,jpds,iwls,iwrs,0,wls2,0)
            iwls = just(lrm,lri)
            iwrs = just(lrm,lrj)
            call prodab(1,0,jpds,iwls,iwrs,0,wls1,0)
          end if
          if (jb_sys > 1) then
            call prodab(1,0,jpdt1,iwlt,iwrt,0,wlt,0)    !bbs_tmp
          end if

        case default
      end select
    end do
    !=========== start  d10_ss(dd,ss,tt)  ==============================
    ! ar(23)-drr(33)-ar(32)  ar(23)-cw(33)-ar(23)
    ! ar(13)-drr(33)-ar(31)  ar(13)-cw(33)-ar(13) w0=1 w1=0
    !if ((lri == 2) .and. (lrj == 4)) then
    !  write(u6,*) 'bbs_tmp err'
    !end if
    wl0 = Zero
    wltmp = Zero
    do lr=lri+1,lrj-1
      list = list3(lri,lrj,lr)
      wltmp = wltmp+(2*vint_ci(list+1)-vint_ci(list))
    end do
    wl0 = wl0-vl0_2*wltmp
    ! drl(33)-bl(23)-ar(32)
    ! drl(33)-bl(13)-ar(31)
    wltmp = Zero
    do lr=1,lri-1
      list = list3(lri,lrj,lr)
      wltmp = wltmp+(vint_ci(list)-2*vint_ci(list+1))
    end do
    wl0 = wl0+vl0_2*wltmp
    ! ar(23)-bl(32)-drl(33)
    ! ar(13)-bl(31)-drl(33)
    wltmp = Zero
    do lr=lrj+1,norb_dz
      list = list3(lri,lrj,lr)
      wltmp = wltmp+(vint_ci(list)-2*vint_ci(list+1))
    end do
    wl0 = wl0+vl0_2*wltmp
    ! ar-arw    arw-al   ar-ar
    list = list3(lri,lrj,lrj)
    wl0 = wl0-vl0_2*vint_ci(list)
    list = list3(lri,lrj,lri)
    wl0 = wl0-vl0_2*(voint(lri,lrj)+vint_ci(list))

    im = Mul(lsm_inn(lri),ns_sm)
    jpdd = 1+im
    iwld = jud(lri)
    iwrd = jud(lrj)
    call prodab(1,0,jpdd,iwld,iwrd,0,wl0,0)
    if (jb_sys > 0) then
      jpdd = jpdd+24
      call prodab(1,0,jpdd,iwld,iwrd,0,wl0,0)
    end if
    do lr=lrj+1,norb_dz
      list = list3(lri,lrj,lr)
      wls = wl0+(vls0-vl0_2)*(vint_ci(list)-2*vint_ci(list+1))-vls1*vint_ci(list)
      ! ar(23)-bl(32)-drl(22) ar(13)-bl(31)-drl(11)
      wlt = wl0+(vlt0-vl0_2)*(vint_ci(list)-2*vint_ci(list+1))-vlt1*vint_ci(list)
      im = Mul(lsm_inn(lri),lsm_inn(lr))
      im = Mul(im,ns_sm)
      jpds = 17+im
      jpdt = 9+im
      iwls = just(lri,lr)
      iwrs = just(lrj,lr)
      iwlt = iwls
      iwrt = iwrs
      call prodab(1,0,jpds,iwls,iwrs,0,wls,0)
      call prodab(1,0,jpdt,iwlt,iwrt,0,wlt,0)
      if (jb_sys > 0) then
        wls_a = wl0+(vls0-vl0_2)*(vint_ci(list)-2*vint_ci(list+1))-vls1_a*vint_ci(list)
        wls_b = -vls10_2b*vint_ci(list)
        iwls = just(lr,lri)
        iwrs = just(lr,lrj)
        call prodab(1,0,jpds,iwls,iwrs,0,wls_a,0)
        iwls = just(lri,lr)
        iwrs = just(lr,lrj)
        call prodab(1,0,jpds,iwls,iwrs,0,wls_b,0)
        iwls = just(lr,lri)
        iwrs = just(lrj,lr)
        call prodab(1,0,jpds,iwls,iwrs,0,wls_b,0)
      end if
      if (jb_sys > 1) then
        jpdt = jpdt+24
        call prodab(1,0,jpdt,iwlt,iwrt,0,wlt,0)
      end if
    end do
    !=== start d5(ss),d40(tt) ==========================================
    do lr=norb_frz+1,lri-1
      list = list3(lri,lrj,lr)
      ! drl(22)-bl(13)-ar(31)
      wls = wl0+(vls0-vl0_2)*(vint_ci(list)-2*vint_ci(list+1))-vls1*vint_ci(list)
      ! drl(22)-bl(23)-ar(32)  drl(11)-bl(13)-ar(31)
      wlt = wl0+(vlt0-vl0_2)*(vint_ci(list)-2*vint_ci(list+1))-vlt1*vint_ci(list)
      im = Mul(lsm_inn(lri),lsm_inn(lr))
      im = Mul(im,ns_sm)
      jpds = 17+im                     !bbs_tmp
      jpdt = 9+im
      iwls = just(lr,lri)
      iwrs = just(lr,lrj)
      iwlt = iwls
      iwrt = iwrs
      call prodab(1,0,jpds,iwls,iwrs,0,wls,0)
      call prodab(1,0,jpdt,iwlt,iwrt,0,wlt,0)
      if (jb_sys > 0) then
        ! drl(11)-bl(23)-ar(32)
        wls = wl0+(vls0-vl0_2)*(vint_ci(list)-2*vint_ci(list+1))-vls1_a*vint_ci(list)
        iwls = just(lri,lr)
        iwrs = just(lrj,lr)
        call prodab(1,0,jpds,iwls,iwrs,0,wls,0)
        ! drl(12)-bl(23)-ar(31) drl(21)-bl(13)-ar(32)
        wls = -vls10_2b*vint_ci(list)
        iwls = just(lri,lr)
        iwrs = just(lr,lrj)
        call prodab(1,0,jpds,iwls,iwrs,0,wls,0)
        iwls = just(lr,lri)
        iwrs = just(lrj,lr)
        call prodab(1,0,jpds,iwls,iwrs,0,wls,0)
      end if
      if (jb_sys > 1) then
        jpat = jpdt+24
        call prodab(1,0,jpat,iwlt,iwrt,0,wlt,0)
      end if
    end do
    !=== end g5,40 =====================================================
  end do
end do

return

end subroutine dbl_space_loop_ijkk_sgezero

subroutine dbl_space_loop_ijkl_sgezero()

use gugaci_global, only: jb_sys, just, lsm_inn, norb_dz, norb_frz, ns_sm, vint_ci
use Symmetry_Info, only: Mul
use Constants, only: Zero, Half
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp) :: im, imi, imij, imik, imil, imj, imk, iml, iwls, iwls1, iwlt, iwrs, iwrs1, iwrt, jpds, jpdt, jpdt1, list, lri, &
                     lrj, lrk, lrl, ni
real(kind=wp) :: db, w0, w1, wls, wls1, wls2, wlt
integer(kind=iwp), external :: list4

!==============================  g11,12  == (v-s)=======================
!==============================  g41,42  == (v-t)=======================
wls1 = Zero
wls2 = Zero
db = jb_sys
do lrl=norb_frz+1,norb_dz-3
  iml = lsm_inn(lrl)
  do lrk=lrl+1,norb_dz-2
    imk = lsm_inn(lrk)
    do lrj=lrk+1,norb_dz-1
      imj = lsm_inn(lrj)
      do lri=norb_dz,lrj+1,-1
        imi = lsm_inn(lri)
        !list = list4(lri,lrj,lrk,lrl)
        list = list4(lrl,lrk,lrj,lri)
        ni = mod(lrk-lrl+lri-lrj,2)
        imik = Mul(imi,imk)
        if (imik == Mul(imj,iml)) then
          im = Mul(imik,ns_sm)
          jpds = 17+im
          jpdt = 9+im
          jpdt1 = jpdt+24
          iwls = just(lrl,lrj)
          iwrs = just(lrk,lri)
          iwlt = iwls
          iwrt = iwrs
          if (jb_sys == 0) then
            ! w0g11=-1/2 w1g11=-3/2                 === g11 ===
            wls = vint_ci(list+1)+vint_ci(list+2)
            ! w0g42=-1/2 w1g42=1/2                  === g42 ===
            wlt = vint_ci(list+1)-vint_ci(list+2)
          end if
          if (jb_sys > 0) then
            ! wog11=-1/2 w1g11=-(db+3)/(2*db+2)
            ! ar(23)bl(32)bl(13)ar(31)
            w1 = -(db+3)/(2*db+2)
            wls = vint_ci(list+1)+(-Half-w1)*vint_ci(list+2)
            ! wog42=-1/2 w1g42=1/2     ar(23)bl(32)bl(23)ar(32)
            wlt = vint_ci(list+1)-vint_ci(list+2)
            ! w0=0 w1=-sqrt(db*(db+2)/(db+1))
            ! ar(13)bl(32)bl(23)ar(31)     ar(23)bl(31)bl(13)ar(32)
            w1 = -sqrt(db*(db+2))/(db+1)
            wls1 = -w1*vint_ci(list+2)
            ! ar(13)bl(31)bl(23)ar(32)      w0=-1/2
            w1 = -(db-1)/(2*db+2)
            wls2 = vint_ci(list+1)+(-Half-w1)*vint_ci(list+2)
            ! ar(13)bl(31)bl(13)ar(31)
            ! w0=-1/2 w1=1/2
            ! wlt1 = wlt
          end if
          !if (ni == 1) wls = -wls
          if (ni == 1) then
            wls = -wls
            wlt = -wlt
            wls1 = -wls1
            wls2 = -wls2
          end if
          call prodab(1,0,jpds,iwls,iwrs,0,wls,0)
          call prodab(1,0,jpdt,iwlt,iwrt,0,wlt,0)
          if (jb_sys > 0) then
            iwrs1 = just(lrk,lri)
            iwls1 = just(lrj,lrl)
            call prodab(1,0,jpds,iwls1,iwrs1,0,wls1,0)
            iwrs1 = just(lri,lrk)
            iwls1 = just(lrl,lrj)
            call prodab(1,0,jpds,iwls1,iwrs1,0,wls1,0)
            iwrs1 = just(lri,lrk)
            iwls1 = just(lrj,lrl)
            call prodab(1,0,jpds,iwls1,iwrs1,0,wls2,0)
          end if
          if (jb_sys > 1) then
            call prodab(1,0,jpdt1,iwlt,iwrt,0,wlt,0)
          end if
        end if
        imil = Mul(imi,iml)
        if (imil == Mul(imj,imk)) then
          im = Mul(imil,ns_sm)
          jpds = 17+im
          jpdt = 9+im
          jpdt1 = jpdt+24
          iwls = just(lrl,lri)
          iwrs = just(lrk,lrj)
          iwlt = iwls
          iwrt = iwrs
          if (jb_sys == 0) then
            ! w0g11=-1/2 w1g11=-3/2                  === g11 ===
            wls = vint_ci(list)+vint_ci(list+1)
            ! w0g42=-1/2 w1g42=1/2                   === g42 ===
            wlt = vint_ci(list+1)-vint_ci(list)
          end if
          if (jb_sys > 0) then
            ! wog11=-1/2 w1g11=-(db+3)/(2*db+2)
            ! ar(23)bl(32)br(31)al(13)
            w1 = -(db+3)/(2*db+2)
            wls = (-Half-w1)*vint_ci(list)+vint_ci(list+1)
            ! wog42=-1/2 w1g42=1/2
            ! ar(23)bl(32)br(23)al(32)
            wlt = vint_ci(list+1)-vint_ci(list)
            ! w0=0 w1=-sqrt(db*(db+2)/(db+1))
            ! ar(23)bl(31)br(32)al(13)
            w1 = -sqrt(db*(db+2))/(db+1)
            wls1 = -w1*vint_ci(list)
            ! ar(13)bl(31)br(32)al(23)
            ! w0=-1/2 w1=-(db-1)/(2*db+2)
            w1 = -(db-1)/(2*db+2)
            wls2 = (-Half-w1)*vint_ci(list)+vint_ci(list+1)
          end if
          if (ni == 1) then
            wls = -wls
            wlt = -wlt
            wls1 = -wls1
            wls2 = -wls2
          end if
          call prodab(1,0,jpds,iwls,iwrs,0,wls,0)
          call prodab(1,0,jpdt,iwlt,iwrt,0,wlt,0)
          if (jb_sys > 0) then
            iwls = just(lri,lrl)
            iwrs = just(lrk,lrj)
            call prodab(1,0,jpds,iwls,iwrs,0,wls1,0)
            iwls = just(lrl,lri)
            iwrs = just(lrj,lrk)
            call prodab(1,0,jpds,iwls,iwrs,0,wls1,0)
            iwls = just(lri,lrl)
            call prodab(1,0,jpds,iwls,iwrs,0,wls2,0)
          end if
          if (jb_sys > 1) then
            call prodab(1,0,jpdt1,iwlt,iwrt,0,wlt,0)
          end if
        end if
        imij = Mul(imi,imj)
        if (imij == Mul(imk,iml)) then
          im = Mul(imij,ns_sm)
          jpds = 17+im
          jpdt = 9+im
          jpdt1 = 24+jpdt
          iwls = just(lrl,lrk)
          iwrs = just(lrj,lri)
          iwlt = iwls
          iwrt = iwrs
          if (jb_sys == 0) then
            ! w0g12=1, w1g12=0                         === g12 ===
            wls = vint_ci(list)+vint_ci(list+2)
            ! w0g41=0, w1g41=1                         === g41 ===
            wlt = vint_ci(list)-vint_ci(list+2)
          end if
          if (jb_sys > 0) then
            w0 = (db+2)/(2*db+2)
            w1 = db/(2*db+2)
            wls = (w0+w1)*vint_ci(list)+(w0-w1)*vint_ci(list+2)
            wlt = vint_ci(list)-vint_ci(list+2)
            w0 = sqrt(db*(db+2))/(2*db+2)
            w1 = -w0
            wls1 = (w0+w1)*vint_ci(list)+(w0-w1)*vint_ci(list+2)
            w0 = db/(2*db+2)
            w1 = (db+2)/(2*db+2)
            wls2 = (w0+w1)*vint_ci(list)+(w0-w1)*vint_ci(list+2)
          end if
          if (ni == 1) then
            wls = -wls
            wlt = -wlt
            wls1 = -wls1
            wls2 = -wls2
          end if
          call prodab(1,0,jpds,iwls,iwrs,0,wls,0)
          call prodab(1,0,jpdt,iwlt,iwrt,0,wlt,0)
          if (jb_sys > 0) then
            iwls = just(lrl,lrk)
            iwrs = just(lri,lrj)
            call prodab(1,0,jpds,iwls,iwrs,0,wls1,0)
            iwls = just(lrk,lrl)
            iwrs = just(lrj,lri)
            call prodab(1,0,jpds,iwls,iwrs,0,wls1,0)
            ! iwls=
            iwrs = just(lri,lrj)
            call prodab(1,0,jpds,iwls,iwrs,0,wls2,0)
          end if
          if (jb_sys > 1) then
            call prodab(1,0,jpdt1,iwlt,iwrt,0,wlt,0)
          end if
        end if
      end do
    end do
  end do
end do

return

end subroutine dbl_space_loop_ijkl_sgezero
