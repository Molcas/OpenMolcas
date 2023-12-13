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

! look for partial loops in active space drt and save them into disk
subroutine guga_ploop(npl,maxplcon)

use gugaci_global, only: idisk_array, idisk_lp, LuLoop, mhlpmax
use Definitions, only: iwp, u6

implicit none
integer(kind=iwp), intent(out) :: npl, maxplcon

maxplcon = 0
mhlpmax = 0
idisk_lp = 0
idisk_array = 0
call idafile(luloop,1,idisk_array,13,idisk_lp)
npl = 0
call vd_lp_search(npl)
call dv_lp_search(npl)
call dd_lp_search(npl)
call dt_lp_search(npl)
call ds_lp_search(npl)
call tv_lp_search(npl)
call td_lp_search(npl)
call tt_lp_search(npl)
call ts_lp_search(npl)
call sv_lp_search(npl)
call sd_lp_search(npl)
call st_lp_search(npl)
call ss_lp_search(npl)

maxplcon = mhlpmax
idisk_lp = 0
call idafile(luloop,1,idisk_array,13,idisk_lp)
write(u6,'(4x,a,i20)') 'Max number of partial loops: ',mhlpmax
!write(u6,*) 'idisk_array'
!write(u6,'(13i8)') idisk_array

return

end subroutine guga_ploop

subroutine sv_lp_search(npl)

use gugaci_global, only: idisk_array, idisk_lp, iml, imr, lp_count, lpblock, lpblock_sv, mhsum, ng_sm
use Definitions, only: iwp, u6

implicit none
integer(kind=iwp), intent(inout) :: npl
integer(kind=iwp) :: iml_

mhsum = 0
lpblock = 0
lp_count(1:22) = 0
idisk_array(10) = idisk_lp

imr = 1
do iml_=1,ng_sm
  iml = iml_
  call act_lp_search(4,4,1)     !id=4
end do
lpblock_sv = lpblock
npl = npl+mhsum
write(u6,'(a15,2i10)') 'sv:',lpblock,mhsum !,idisk_array(10),idisk_l

return

end subroutine sv_lp_search

subroutine sd_lp_search(npl)

use gugaci_global, only: idisk_array, idisk_lp, iml, imr, lp_count, lpblock, lpblock_sd, mhsum, ng_sm
use Definitions, only: iwp, u6

implicit none
integer(kind=iwp), intent(inout) :: npl
integer(kind=iwp) :: iml_, imr_

mhsum = 0
lpblock = 0
lp_count(1:22) = 0
idisk_array(11) = idisk_lp

do iml_=1,ng_sm
  iml = iml_
  do imr_=1,ng_sm
    imr = imr_
    call act_lp_search(1,4,2)    !id=1
  end do
end do
lpblock_sd = lpblock
npl = npl+mhsum
write(u6,'(a15,2i10)') 'sd:',lpblock,mhsum

return

end subroutine sd_lp_search

subroutine st_lp_search(npl)

use gugaci_global, only: idisk_array, idisk_lp, iml, imr, lp_count, lpblock, lpblock_st, mhsum, ng_sm
use Definitions, only: iwp, u6

implicit none
integer(kind=iwp), intent(inout) :: npl
integer(kind=iwp) :: iml_, imr_

mhsum = 0
lpblock = 0
lp_count(1:22) = 0
idisk_array(12) = idisk_lp

do iml_=1,ng_sm
  iml = iml_
  do imr_=1,ng_sm
    imr = imr_
    call act_lp_search(3,4,3)   !id=3
  end do
end do
lpblock_st = lpblock
npl = npl+mhsum
write(u6,'(a15,2i10)') 'st:',lpblock,mhsum

return

end subroutine st_lp_search

subroutine ss_lp_search(npl)

use gugaci_global, only: idisk_array, idisk_lp, iml, imr, lp_count, lpblock, lpblock_ss, mhsum, ng_sm
use Definitions, only: iwp, u6

implicit none
integer(kind=iwp), intent(inout) :: npl
integer(kind=iwp) :: iml_, imr_

mhsum = 0
lpblock = 0
lp_count(1:22) = 0
idisk_array(13) = idisk_lp

do iml_=1,ng_sm
  iml = iml_
  do imr_=1,ng_sm
    imr = imr_
    call act_lp_search(3,4,4)    !id=3
  end do
end do
lpblock_ss = lpblock
npl = npl+mhsum
write(u6,'(a15,2i10)') 'ss:',lpblock,mhsum

return

end subroutine ss_lp_search

subroutine tv_lp_search(npl)

use gugaci_global, only: idisk_array, idisk_lp, iml, imr, lp_count, lpblock, lpblock_tv, mhsum, ng_sm
use Definitions, only: iwp, u6

implicit none
integer(kind=iwp), intent(inout) :: npl
integer(kind=iwp) :: iml_

mhsum = 0
lpblock = 0
lp_count(1:22) = 0
idisk_array(6) = idisk_lp

imr = 1
do iml_=1,ng_sm
  iml = iml_
  call act_lp_search(4,3,1)       !id=4
end do
lpblock_tv = lpblock
npl = npl+mhsum
write(u6,'(a15,2i10)') 'tv:',lpblock,mhsum

return

end subroutine tv_lp_search

subroutine td_lp_search(npl)

use gugaci_global, only: idisk_array, idisk_lp, iml, imr, lp_count, lpblock, lpblock_td, mhsum, ng_sm
use Definitions, only: iwp, u6

implicit none
integer(kind=iwp), intent(inout) :: npl
integer(kind=iwp) :: iml_, imr_

mhsum = 0
lpblock = 0
lp_count(1:22) = 0
idisk_array(7) = idisk_lp

do iml_=1,ng_sm
  iml = iml_
  do imr_=1,ng_sm
    imr = imr_
    call act_lp_search(1,3,2)    !id=1
  end do
end do
lpblock_td = lpblock
npl = npl+mhsum
write(u6,'(a15,2i10)') 'td:',lpblock,mhsum

return

end subroutine td_lp_search

subroutine tt_lp_search(npl)

use gugaci_global, only: idisk_array, idisk_lp, iml, imr, lp_count, lpblock, lpblock_tt, mhsum, ng_sm
use Definitions, only: iwp, u6

implicit none
integer(kind=iwp), intent(inout) :: npl
integer(kind=iwp) :: iml_, imr_

mhsum = 0
lpblock = 0
lp_count(1:22) = 0
idisk_array(8) = idisk_lp

do iml_=1,ng_sm
  iml = iml_
  do imr_=1,ng_sm
    imr = imr_
    call act_lp_search(3,3,3)    !id=3
  end do
end do
lpblock_tt = lpblock
npl = npl+mhsum
write(u6,'(a15,2i10)') 'tt:',lpblock,mhsum

return

end subroutine tt_lp_search

subroutine ts_lp_search(npl)

use gugaci_global, only: idisk_array, idisk_lp, iml, imr, lp_count, lpblock, lpblock_ts, mhsum, ng_sm
use Definitions, only: iwp, u6

implicit none
integer(kind=iwp), intent(inout) :: npl
integer(kind=iwp) :: iml_, imr_

mhsum = 0
lpblock = 0
lp_count(1:22) = 0
idisk_array(9) = idisk_lp

do iml_=1,ng_sm
  iml = iml_
  do imr_=1,ng_sm
    imr = imr_
    call act_lp_search(3,3,4)    !id=3
  end do
end do
lpblock_ts = lpblock
npl = npl+mhsum
write(u6,'(a15,2i10)') 'ts:',lpblock,mhsum

return

end subroutine ts_lp_search

subroutine dv_lp_search(npl)

use gugaci_global, only: idisk_array, idisk_lp, iml, imr, lp_count, lpblock, lpblock_dv, mhsum, ng_sm
use Definitions, only: iwp, u6

implicit none
integer(kind=iwp), intent(inout) :: npl
integer(kind=iwp) :: iml_

mhsum = 0
lpblock = 0
lp_count(1:22) = 0
idisk_array(2) = idisk_lp

imr = 1
do iml_=1,ng_sm
  iml = iml_
  call act_lp_search(1,2,1)      !id=1
end do
lpblock_dv = lpblock
npl = npl+mhsum
write(u6,'(a15,2i10)') 'dv:',lpblock,mhsum

return

end subroutine dv_lp_search

subroutine dd_lp_search(npl)

use gugaci_global, only: idisk_array, idisk_lp, iml, imr, lp_count, lpblock, lpblock_dd, mhsum, ng_sm
use Definitions, only: iwp, u6

implicit none
integer(kind=iwp), intent(inout) :: npl
integer(kind=iwp) :: iml_, imr_

mhsum = 0
lpblock = 0
lp_count(1:22) = 0
idisk_array(3) = idisk_lp

do iml_=1,ng_sm
  iml = iml_
  do imr_=1,ng_sm
    imr = imr_
    call act_lp_search(3,2,2)    !id=3
  end do
end do
lpblock_dd = lpblock
npl = npl+mhsum
write(u6,'(a15,2i10)') 'dd:',lpblock,mhsum

return

end subroutine dd_lp_search

subroutine dt_lp_search(npl)

use gugaci_global, only: idisk_array, idisk_lp, iml, imr, lp_count, lpblock, lpblock_dt, mhsum, ng_sm
use Definitions, only: iwp, u6

implicit none
integer(kind=iwp), intent(inout) :: npl
integer(kind=iwp) :: iml_, imr_

mhsum = 0
lpblock = 0
lp_count(1:22) = 0
idisk_array(4) = idisk_lp

do iml_=1,ng_sm
  iml = iml_
  do imr_=1,ng_sm
    imr = imr_
    call act_lp_search(2,2,3)      !id=2
  end do
end do
lpblock_dt = lpblock
npl = npl+mhsum
write(u6,'(a15,2i10)') 'dt:',lpblock,mhsum

return

end subroutine dt_lp_search

subroutine ds_lp_search(npl)

use gugaci_global, only: idisk_array, idisk_lp, iml, imr, lp_count, lpblock, lpblock_ds, mhsum, ng_sm
use Definitions, only: iwp, u6

implicit none
integer(kind=iwp), intent(inout) :: npl
integer(kind=iwp) :: iml_, imr_

mhsum = 0
lpblock = 0
lp_count(1:22) = 0
idisk_array(5) = idisk_lp

do iml_=1,ng_sm
  iml = iml_
  do imr_=1,ng_sm
    imr = imr_
    call act_lp_search(2,2,4)    !id=2
  end do
end do
lpblock_ds = lpblock
npl = npl+mhsum
write(u6,'(a15,2i10)') 'ds:',lpblock,mhsum

return

end subroutine ds_lp_search

subroutine vd_lp_search(npl)

use gugaci_global, only: idisk_array, idisk_lp, iml, imr, lp_count, lpblock, lpblock_vd, mhsum, ng_sm
use Definitions, only: iwp, u6

implicit none
integer(kind=iwp), intent(inout) :: npl
integer(kind=iwp) :: imr_

mhsum = 0
lpblock = 0
lp_count(1:22) = 0
idisk_array(1) = idisk_lp

iml = 1
do imr_=1,ng_sm
  imr = imr_
  call act_lp_search(2,1,2)       !id=2
end do
lpblock_vd = lpblock
npl = npl+mhsum
write(u6,'(a15,2i10)') 'vd:',lpblock,mhsum

return

end subroutine vd_lp_search

subroutine act_lp_search(id,iptyl,iptyr)

use gugaci_global, only: iml, imr, ipae, ipael, jml, jmr, jpad, jpadl, jpadlr, jpae, jpael, map_jplr, ndim, ng_sm, ns_sm, nu_ad, &
                         nu_ae
use Definitions, only: iwp

implicit none
integer(kind=iwp), intent(in) :: id, iptyl, iptyr
integer(kind=iwp) :: ide, ipaer, jml_, jmlend, jmlsta, jmr_, jmrend, jmrsta, jpadr, jpaer, jptyl, jptyr, jptyrend

call get_jp(iptyl,iml,ipael,0)
call get_jp(iptyr,imr,ipaer,0)

jpael = nu_ae(ipael)
jpaer = nu_ae(ipaer)
if ((jpael == 0) .or. (jpaer == 0)) return
ide = 0
if ((iptyl == 4) .and. (iptyr == 4)) ide = 1
if ((iptyl == 4) .and. (iptyr == 3)) ide = 2
if ((iptyl == 3) .and. (iptyr == 4)) ide = 2
do jptyl=1,6            !1,6???
  jptyrend = 6
  if (jptyl == 1) jptyrend = 1
  jmlsta = 1
  jmlend = ng_sm
  if (jptyl == 1) jmlsta = ns_sm
  if (jptyl == 1) jmlend = ns_sm
  do jml_=jmlsta,jmlend
    jml = jml_          ! jml is in global module, is this necessary?
    call get_jp(jptyl,jml,jpadl,0)
    if (nu_ad(jpadl) == 0) cycle
    jpad = jpadl
    jpae = jpael
    ipae = ipael
    call seg_drt()
    if (ndim == 0) cycle
    call copy_to_drtl()
    do jptyr=1,jptyrend
      jmrsta = 1
      jmrend = ng_sm
      if (jptyr == 1) jmrsta = ns_sm
      if (jptyr == 1) jmrend = ns_sm
      do jmr_=jmrsta,jmrend
        jmr = jmr_      ! jmr is in global module, is this necessary?
        call get_jp(jptyr,jmr,jpadr,0)
        if (nu_ad(jpadr) == 0) cycle
        jpad = jpadr
        jpae = jpaer
        ipae = ipaer
        call seg_drt()
        if (ndim == 0) cycle
        jpadlr = map_jplr(jptyl,jptyr)
        if (id == 1) call lp_head_in_dbl_1()
        if (id == 2) call lp_head_in_dbl_2()
        if (id == 3) call lp_head_in_dbl_3(ide)
        if (id == 4) call lp_head_in_dbl_4()
        if (jpadl /= jpadr) cycle
        if (id == 1) call lp_head_in_act_1()
        if (id == 2) call lp_head_in_act_2()
        if (id == 3) call lp_head_in_act_3(ide)
        if (id == 4) call lp_head_in_act_4()
      end do
    end do
  end do
end do

return

end subroutine act_lp_search

subroutine lp_head_in_act_1()      !for dv,td,sd

use gugaci_global, only: iml, imr, intind_ijka, linelp, log_prod, logic_br, lpblock, lsm_inn, mhlp, ngw2, ngw3, nlg1, nlg2, &
                         norb_dz, norb_frz, norb_inn
use Symmetry_Info, only: Mul
use Definitions, only: iwp

implicit none
integer(kind=iwp) :: ijk, imlr, intpos, lma, lmai, lmk, lra, lrai, lraj, lrak, lsmi, lsmij, lsmj, lsmk, mh

imlr = Mul(iml,imr)
do lra=norb_dz+1,norb_inn
  lma = lsm_inn(lra)
  if (lma /= imlr) cycle
  ! line=1 a&r<-->a^r
  call head_ar_at_given_orb(mh,lra)
  call link_c1_to_given_orb_coe(mh,lra+1,norb_inn)
  if (mh == 0) cycle
  call value_sort_ploop(mh,.true.,.true.,.false.)
  if (log_prod == 3) then
    linelp = 1
    mhlp = mh
    nlg1 = lra
    nlg2 = 0
    call ext_head_in_act()
  else
    call save_lp(1,mh,lra,0)
    lpblock = lpblock+1
  end if
end do

do lrai=norb_dz+1,norb_inn
  lsmi = lsm_inn(lrai)
  do lraj=lrai+1,norb_inn
    lsmj = lsm_inn(lraj)
    lsmij = Mul(lsmi,lsmj)
    lmk = Mul(lsmij,imlr)
    do lrak=lraj+1,norb_inn
      lsmk = lsm_inn(lrak)
      if (lmk /= lsmk) cycle
      ijk = lrai-norb_frz+ngw2(lraj-norb_frz)+ngw3(lrak-norb_frz)
      intpos = intind_ijka(ijk)
      ! line=4  a&r--b&r--b^r<-->a^r
      call head_ar_at_given_orb(mh,lrai)
      call link_c1_to_given_orb(mh,lrai+1,lraj-1)
      call link_b4_at_given_orb(mh)
      logic_br(1:mh) = .true.
      call link_c2_to_given_orb(mh,lraj+1,lrak-1)
      call link_b2_at_given_orb(mh)
      call link_c1_to_given_orb(mh,lrak+1,norb_inn)
      if (mh /= 0) then
        call value_sort_ploop(mh,.false.,.true.,.true.)
        if (log_prod == 3) then
          linelp = 4
          mhlp = mh
          nlg1 = intpos
          nlg2 = 0
          call ext_head_in_act()
        else
          call save_lp(4,mh,intpos,0)
          lpblock = lpblock+1
        end if
      end if
      ! line=7 a&r--b&l--b^l<-->a^r
      call head_ar_at_given_orb(mh,lrai)
      call link_c1_to_given_orb(mh,lrai+1,lraj-1)
      call link_b3_at_given_orb(mh)
      logic_br(1:mh) = .true.
      call link_c2_to_given_orb(mh,lraj+1,lrak-1)
      call link_b1_at_given_orb(mh)
      call link_c1_to_given_orb(mh,lrak+1,norb_inn)
      if (mh == 0) cycle
      call value_sort_ploop(mh,.false.,.true.,.true.)
      if (log_prod == 3) then
        linelp = 7
        mhlp = mh
        nlg1 = intpos
        nlg2 = 0
        call ext_head_in_act()
      else
        call save_lp(7,mh,intpos,0)
        lpblock = lpblock+1
      end if
    end do
  end do
end do
do lrai=norb_dz+2,norb_inn
  lmai = lsm_inn(lrai)
  if (lmai /= imlr) cycle
  do lrak=norb_dz+1,lrai-1
    ! line=10 d&rr--b^r<-->a^r
    call head_drr_at_given_orb(mh,lrak)
    logic_br(1:mh) = .true.
    call link_c2_to_given_orb(mh,lrak+1,lrai-1)
    call link_b2_at_given_orb(mh)
    call link_c1_to_given_orb(mh,lrai+1,norb_inn)
    if (mh /= 0) then
      call value_sort_ploop(mh,.false.,.true.,.true.)
      if (log_prod == 3) then
        linelp = 10
        mhlp = mh
        nlg1 = lrak
        nlg2 = lrai
        call ext_head_in_act()
      else
        call save_lp(10,mh,lrak,lrai)
        lpblock = lpblock+1
      end if
    end if
    ! line=12 d&rl--b^l<-->a^r
    call head_drl_at_given_orb(mh,lrak)
    call link_c2_to_given_orb(mh,lrak+1,lrai-1)
    call link_b1_at_given_orb(mh)
    call link_c1_to_given_orb(mh,lrai+1,norb_inn)
    if (mh == 0) cycle
    call value_sort_ploop(mh,.false.,.true.,.true.)
    if (log_prod == 3) then
      linelp = 12
      mhlp = mh
      nlg1 = lrak
      nlg2 = lrai
      call ext_head_in_act()
    else
      call save_lp(12,mh,lrak,lrai)
      lpblock = lpblock+1
    end if
  end do
end do

return

end subroutine lp_head_in_act_1

subroutine lp_head_in_act_2()      !for vd,dt,ds

use gugaci_global, only: iml, imr, intind_ijka, linelp, log_prod, logic_br, lpblock, lsm_inn, mhlp, ngw2, ngw3, nlg1, nlg2, &
                         norb_dz, norb_frz, norb_inn
use Symmetry_Info, only: Mul
use Definitions, only: iwp

implicit none
integer(kind=iwp) :: ijk, imlr, intpos, lma, lmai, lmk, lra, lrai, lraj, lrak, lrd, lsmi, lsmij, lsmj, lsmk, mh

imlr = Mul(iml,imr)

do lra=norb_dz+1,norb_inn
  lma = lsm_inn(lra)
  if (lma /= imlr) cycle
  do lrd=lra+1,norb_inn
    ! line=2 a&r--d&l^r<-->a^l
    call head_ar_at_given_orb(mh,lra)
    call link_c1_to_given_orb(mh,lra+1,lrd-1)
    call link_d10_at_given_orb(mh)
    call link_c1_to_given_orb(mh,lrd+1,norb_inn)
    if (mh == 0) cycle
    call value_sort_ploop(mh,.false.,.true.,.false.)
    if (log_prod == 3) then
      linelp = 2
      mhlp = mh
      nlg1 = lra
      nlg2 = lrd
      !call print_lp()
      call ext_head_in_act()
    else
      linelp = 2
      mhlp = mh
      nlg1 = lra
      nlg2 = lrd
      !call print_lp()
      call save_lp(2,mh,lra,lrd)
      lpblock = lpblock+1
    end if
  end do
end do

do lrai=norb_dz+1,norb_inn
  lsmi = lsm_inn(lrai)
  do lraj=lrai+1,norb_inn
    lsmj = lsm_inn(lraj)
    lsmij = Mul(lsmi,lsmj)
    lmk = Mul(lsmij,imlr)
    do lrak=lraj+1,norb_inn
      lsmk = lsm_inn(lrak)
      if (lmk /= lsmk) cycle
      ijk = lrai-norb_frz+ngw2(lraj-norb_frz)+ngw3(lrak-norb_frz)
      intpos = intind_ijka(ijk)
      ! line=6  a&r--b&l--b^r<-->a^l
      call head_ar_at_given_orb(mh,lrai)
      call link_c1_to_given_orb(mh,lrai+1,lraj-1)
      call link_b3_at_given_orb(mh)
      logic_br(1:mh) = .true.
      call link_c2_to_given_orb(mh,lraj+1,lrak-1)
      call link_b2_at_given_orb(mh)
      call link_c1_to_given_orb(mh,lrak+1,norb_inn)
      if (mh == 0) cycle
      call value_sort_ploop(mh,.false.,.true.,.true.)
      if (log_prod == 3) then
        linelp = 6
        mhlp = mh
        nlg1 = intpos
        nlg2 = 0
        !call print_lp()
        call ext_head_in_act()
      else
        call save_lp(6,mh,intpos,0)
        lpblock = lpblock+1
      end if
    end do
  end do
end do

do lrai=norb_dz+2,norb_inn
  lmai = lsm_inn(lrai)
  if (lmai /= imlr) cycle
  do lrak=norb_dz+1,lrai-1
    ! line=11 d&rl--b^r<-->a^l
    call head_drl_at_given_orb(mh,lrak)
    call link_c2_to_given_orb(mh,lrak+1,lrai-1)
    call link_b2_at_given_orb(mh)
    call link_c1_to_given_orb(mh,lrai+1,norb_inn)
    if (mh == 0) cycle
    call value_sort_ploop(mh,.false.,.true.,.true.)
    if (log_prod == 3) then
      linelp = 11
      mhlp = mh
      nlg1 = lrak
      nlg2 = lrai
      !call print_lp()
      call ext_head_in_act()
    else
      call save_lp(11,mh,lrak,lrai)
      lpblock = lpblock+1
    end if
  end do
end do

return

end subroutine lp_head_in_act_2

subroutine lp_head_in_act_3(ide)      !for ide=0:dd,tt,ide=1:ss,id

use gugaci_global, only: iml, imr, logic_br, lpblock, lsm_inn, norb_dz, norb_inn
use Symmetry_Info, only: Mul
use Definitions, only: iwp

implicit none
integer(kind=iwp), intent(in) :: ide
integer(kind=iwp) :: imlr, lra, lrai, lraj, lsmi, lsmij, lsmj, mh

imlr = Mul(iml,imr)

do lra=norb_dz+1,norb_inn
  ! line=9 d&r&l-
  call head_drl_at_given_orb(mh,lra)
  call link_c2_to_given_orb(mh,lra+1,norb_inn)
  if (mh == 0) cycle
  call value_sort_ploop(mh,.false.,.true.,.true.)
  call save_lp(9,mh,lra,0)
  lpblock = lpblock+1
end do
do lrai=norb_dz+1,norb_inn
  lsmi = lsm_inn(lrai)
  do lraj=lrai+1,norb_inn
    lsmj = lsm_inn(lraj)
    lsmij = Mul(lsmi,lsmj)
    if (lsmij /= imlr) cycle
    ! line=5 a&r-b&l-
    call head_ar_at_given_orb(mh,lrai)
    call link_c1_to_given_orb(mh,lrai+1,lraj-1)
    call link_b3_at_given_orb(mh)
    logic_br(1:mh) = .true.
    call link_c2_to_given_orb(mh,lraj+1,norb_inn)
    if (mh == 0) cycle
    select case (ide)
      case (0)
        call value_sort_ploop(mh,.false.,.true.,.true.)
      case default ! (1)
        call value_sort_ploop(mh,.false.,.true.,.false.)    !ss
      case (2)
        call value_sort_ploop(mh,.false.,.false.,.true.)    !st,ts
    end select
    call save_lp(5,mh,lrai,lraj)
    lpblock = lpblock+1
  end do
end do

return

end subroutine lp_head_in_act_3

subroutine lp_head_in_act_4()

use gugaci_global, only: iml, imr, linelp, log_prod, logic_br, lpblock, lsm_inn, mhlp, nlg1, nlg2, norb_dz, norb_inn
use Symmetry_Info, only: Mul
use Definitions, only: iwp

implicit none
integer(kind=iwp) :: imlr, lra, lrai, lraj, lsmi, lsmij, lsmj, mh

imlr = Mul(iml,imr)

do lra=norb_dz+1,norb_inn
  ! line=8 d&rr-
  call head_drr_at_given_orb(mh,lra)
  logic_br(1:mh) = .true.
  call link_c2_to_given_orb(mh,lra+1,norb_inn)
  if (mh == 0) cycle
  call value_sort_ploop(mh,.false.,.true.,.true.)
  if (log_prod == 3) then
    linelp = 8
    mhlp = mh
    nlg1 = lra
    nlg2 = 0
    call ext_head_in_act()
  else
    call save_lp(8,mh,lra,0)
    lpblock = lpblock+1
  end if
end do
do lrai=norb_dz+1,norb_inn
  lsmi = lsm_inn(lrai)
  do lraj=lrai+1,norb_inn
    lsmj = lsm_inn(lraj)
    lsmij = Mul(lsmi,lsmj)
    if (lsmij /= imlr) cycle
    call head_ar_at_given_orb(mh,lrai)
    call link_c1_to_given_orb(mh,lrai+1,lraj-1)
    call link_b4_at_given_orb(mh)
    logic_br(1:mh) = .true.
    call link_c2_to_given_orb(mh,lraj+1,norb_inn)
    if (mh == 0) cycle
    call value_sort_ploop(mh,.false.,.true.,.true.)
    if (log_prod == 3) then
      linelp = 3
      mhlp = mh
      nlg1 = lrai
      nlg2 = lraj
      call ext_head_in_act()
    else
      call save_lp(3,mh,lrai,lraj)
      lpblock = lpblock+1
    end if
  end do
end do

return

end subroutine lp_head_in_act_4

subroutine lp_head_in_dbl_1()      !for dv,sd,td

use gugaci_global, only: iml, imr, jml, jmr, logic_br, lpblock, lsm_inn, ngw2, ngw3, norb_dz, norb_frz, norb_inn
use Symmetry_Info, only: Mul
use Definitions, only: iwp

implicit none
integer(kind=iwp) :: imlr, jk, jmlr, lend, lra, lrai, lraj, lsma, lsmact, lsmi, lsmij, lsmj, lsta, mh

imlr = Mul(iml,imr)
jmlr = Mul(jml,jmr)
lsmact = Mul(imlr,jmlr)
lsta = norb_dz+1
lend = norb_inn

if (imlr == jmlr) then       !!! imlr == 1

  ! line=13 -c'-
  call link_c1_to_given_orb_coe(mh,lsta,lend)
  if (mh /= 0) then
    call value_sort_ploop(mh,.true.,.true.,.true.)
    call save_lp(13,mh,0,0)
    lpblock = lpblock+1
  end if
end if

! line=15 -b^l-
do lra=norb_dz+1,norb_inn
  lsma = lsm_inn(lra)

  if (lsma /= lsmact) cycle

  if (jml == jmr) then
    logic_br(1) = .false.
    call link_c2_to_given_orb(mh,norb_dz+1,lra-1)
    call link_b1_at_given_orb(mh)
    call link_c1_to_given_orb(mh,lra+1,norb_inn)
    if (mh /= 0) then
      call value_sort_ploop(mh,.false.,.true.,.true.)
      call save_lp(15,mh,lra,1)
      lpblock = lpblock+1
    end if
  end if
  logic_br(1) = .true.
  call link_c2_to_given_orb(mh,norb_dz+1,lra-1)
  call link_b1_at_given_orb(mh)        !b^l
  call link_c1_to_given_orb(mh,lra+1,norb_inn)
  if (mh /= 0) then
    call value_sort_ploop(mh,.false.,.true.,.true.)
    call save_lp(15,mh,lra,2)
    lpblock = lpblock+1
  end if

  ! line=16 -b^r-
  if (jml == jmr) then
    logic_br(1) = .false.
    call link_c2_to_given_orb(mh,norb_dz+1,lra-1)
    call link_b2_at_given_orb(mh)
    call link_c1_to_given_orb(mh,lra+1,norb_inn)
    if (mh /= 0) then
      call value_sort_ploop(mh,.false.,.true.,.true.)
      call save_lp(16,mh,lra,1)
      lpblock = lpblock+1
    end if
  end if

  logic_br(1) = .true.
  call link_c2_to_given_orb(mh,norb_dz+1,lra-1)
  call link_b2_at_given_orb(mh)
  call link_c1_to_given_orb(mh,lra+1,norb_inn)
  if (mh == 0) cycle
  call value_sort_ploop(mh,.false.,.true.,.true.)
  call save_lp(16,mh,lra,2)
  lpblock = lpblock+1

end do

! line=20 -b&r-b^r-
do lrai=norb_dz+1,norb_inn-1
  lsmi = lsm_inn(lrai)
  do lraj=lrai+1,norb_inn
    lsmj = lsm_inn(lraj)
    lsmij = Mul(lsmi,lsmj)
    if (lsmij /= lsmact) cycle
    jk = ngw2(lrai-norb_frz)+ngw3(lraj-norb_frz)

    call link_c1_to_given_orb(mh,norb_dz+1,lrai-1)
    call link_b4_at_given_orb(mh)        !b&r
    logic_br(1:mh) = .true.
    call link_c2_to_given_orb(mh,lrai+1,lraj-1)
    call link_b2_at_given_orb(mh)        !b^r
    call link_c1_to_given_orb(mh,lraj+1,norb_inn)
    if (mh /= 0) then
      call value_sort_ploop(mh,.false.,.true.,.true.)
      call save_lp(20,mh,jk,0)
      lpblock = lpblock+1
    end if
    ! line=22 -b&l-b^l-
    call link_c1_to_given_orb(mh,norb_dz+1,lrai-1)
    call link_b3_at_given_orb(mh)        !b&l
    logic_br(1:mh) = .true.
    call link_c2_to_given_orb(mh,lrai+1,lraj-1)
    call link_b1_at_given_orb(mh)        !b^l
    call link_c1_to_given_orb(mh,lraj+1,norb_inn)
    if (mh == 0) cycle
    call value_sort_ploop(mh,.false.,.true.,.true.)
    call save_lp(22,mh,jk,0)
    lpblock = lpblock+1
  end do
end do

return

end subroutine lp_head_in_dbl_1

!subroutine lp_head_in_dbl_1_mrpt2()    !for dv,sd,td       !200709
!
!use gugaci_global, only: iml, imr, jml, jmr, jpadlr, jpadlrel, linelp, logic_br, lsm_inn, mhlp, ngw2, ngw3, nlg1, nlg2, norb_dz, &
!                         norb_frz, norb_inn
!use Symmetry_Info, only: Mul
!use Definitions, only: iwp
!
!implicit none
!integer(kind=iwp) :: imlr, jk, jmlr, lend, lra, lrai, lraj, lsma, lsmact, lsmi, lsmij, lsmj, lsta, mh
!logical(kind=iwp) :: do_15, do_16
!
!imlr = Mul(iml,imr)
!jmlr = Mul(jml,jmr)
!lsmact = Mul(imlr,jmlr)
!lsta = norb_dz+1
!lend = norb_inn
!
!do_15 = .false.
!do_16 = .false.
!select case (jpadlrel(jpadlr)) !200709
!  case default ! (1)
!    if (imlr == jmlr) then               !!! imlr == 1
!      ! line=13 -c'-
!      call link_c1_to_given_orb_coe(mh,lsta,lend)
!      if (mh == 0) then
!        do_15 = .true.
!      else
!        call value_sort_ploop(mh,.true.,.true.,.true.)
!        linelp = 13
!        mhlp = mh
!        nlg1 = 0
!        nlg2 = 0
!        call ext_head_in_dbl()
!      end if
!    end if
!
!    if (.not. do_15) then
!      ! line=20 -b&r-b^r-
!      do lrai=norb_dz+1,norb_inn-1
!        lsmi = lsm_inn(lrai)
!        do lraj=lrai+1,norb_inn
!          lsmj = lsm_inn(lraj)
!          lsmij = Mul(lsmi,lsmj)
!          if (lsmij /= lsmact) cycle
!          jk = ngw2(lrai-norb_frz)+ngw3(lraj-norb_frz)
!
!          call link_c1_to_given_orb(mh,norb_dz+1,lrai-1)
!          call link_b4_at_given_orb(mh)        !b&r
!          logic_br(1:mh) = .true.
!          call link_c2_to_given_orb(mh,lrai+1,lraj-1)
!          call link_b2_at_given_orb(mh)        !b^r
!          call link_c1_to_given_orb(mh,lraj+1,norb_inn)
!          if (mh /= 0) then
!            call value_sort_ploop(mh,.false.,.true.,.true.)
!            linelp = 20
!            mhlp = mh
!            nlg1 = jk
!            nlg2 = 0
!            call ext_head_in_dbl()
!          end if
!          ! line=22 -b&l-b^l-
!          call link_c1_to_given_orb(mh,norb_dz+1,lrai-1)
!          call link_b3_at_given_orb(mh)        !b&l
!          logic_br(1:mh) = .true.
!          call link_c2_to_given_orb(mh,lrai+1,lraj-1)
!          call link_b1_at_given_orb(mh)        !b^l
!          call link_c1_to_given_orb(mh,lraj+1,norb_inn)
!          if (mh == 0) cycle
!          call value_sort_ploop(mh,.false.,.true.,.true.)
!          linelp = 22
!          mhlp = mh
!          nlg1 = jk
!          nlg2 = 0
!          call ext_head_in_dbl()
!        end do
!      end do
!    end if
!
!  case (2)
!
!  case (3)
!    do_15 = .true.
!
!  case (4)
!    do_16 = .true.
!end select
!
!! line=15 -b^l-
!if (do_15) then
!  do lra=norb_dz+1,norb_inn
!    lsma = lsm_inn(lra)
!
!    if (lsma /= lsmact) cycle
!
!    if (jml == jmr) then
!      logic_br(1) = .false.
!      call link_c2_to_given_orb(mh,norb_dz+1,lra-1)
!      call link_b1_at_given_orb(mh)
!      call link_c1_to_given_orb(mh,lra+1,norb_inn)
!      if (mh /= 0) then
!        call value_sort_ploop(mh,.false.,.true.,.true.)
!        linelp = 15
!        mhlp = mh
!        nlg1 = lra
!        nlg2 = 1
!        call ext_head_in_dbl()
!        !write(u6,*) mh,linelp
!      end if
!    end if
!    logic_br(1) = .true.
!    call link_c2_to_given_orb(mh,norb_dz+1,lra-1)
!    call link_b1_at_given_orb(mh)        !b^l
!    call link_c1_to_given_orb(mh,lra+1,norb_inn)
!    if (mh == 0) then
!      do_16 = .true.
!      exit
!    else
!      call value_sort_ploop(mh,.false.,.true.,.true.)
!      linelp = 15
!      mhlp = mh
!      nlg1 = lra
!      nlg2 = 2
!      call ext_head_in_dbl()
!      !write(u6,*) mh,linelp
!    end if
!  end do
!end if
!
!! line=16 -b^r-
!if (do_16) then
!  do lra=norb_dz+1,norb_inn
!    lsma = lsm_inn(lra)
!    if (jml == jmr) then
!      logic_br(1) = .false.
!      call link_c2_to_given_orb(mh,norb_dz+1,lra-1)
!      call link_b2_at_given_orb(mh)
!      call link_c1_to_given_orb(mh,lra+1,norb_inn)
!      if (mh /= 0) then
!        call value_sort_ploop(mh,.false.,.true.,.true.)
!        linelp = 16
!        mhlp = mh
!        nlg1 = lra
!        nlg2 = 1
!        call ext_head_in_dbl()
!        !write(u6,*) mh,linelp
!      end if
!    end if
!
!    logic_br(1) = .true.
!    call link_c2_to_given_orb(mh,norb_dz+1,lra-1)
!    call link_b2_at_given_orb(mh)
!    call link_c1_to_given_orb(mh,lra+1,norb_inn)
!    if (mh == 0) cycle
!    call value_sort_ploop(mh,.false.,.true.,.true.)
!    linelp = 16
!    mhlp = mh
!    nlg1 = lra
!    nlg2 = 2
!    call ext_head_in_dbl()
!    !write(u6,*) mh,linelp
!
!  end do
!end if
!
!return
!
!end subroutine lp_head_in_dbl_1_mrpt2

subroutine lp_head_in_dbl_2()      !for vd,dt,ds

use gugaci_global, only: iml, imr, jml, jmr, linelp, log_prod, logic_br, lpblock, lsm_inn, mhlp, ngw2, ngw3, nlg1, nlg2, norb_dz, &
                         norb_frz, norb_inn
use Symmetry_Info, only: Mul
use Definitions, only: iwp

implicit none
integer(kind=iwp) :: imlr, jk, jmlr, lend, lra, lrai, lraj, lsma, lsmact, lsmi, lsmij, lsmj, lsta, mh

imlr = Mul(iml,imr)
jmlr = Mul(jml,jmr)
lsmact = Mul(imlr,jmlr)
lsta = norb_dz+1
lend = norb_inn

if (imlr == jmlr) then       !!! imlr == 1

  ! line=13 -c'-
  call link_c1_to_given_orb_coe(mh,lsta,lend)
  if (mh /= 0) then
    call value_sort_ploop(mh,.true.,.true.,.true.)
    if (log_prod == 3) then
      ! mrpt2, do not save loop
      linelp = 13
      mhlp = mh
      nlg1 = 0
      nlg2 = 0
      call ext_head_in_dbl()
    else
      call save_lp(13,mh,0,0)
      lpblock = lpblock+1
    end if
  end if
end if

! line=16 -b^r-
do lra=norb_dz+1,norb_inn
  lsma = lsm_inn(lra)

  if (lsma == lsmact) then

    if (jml == jmr) then
      logic_br(1) = .false.
      call link_c2_to_given_orb(mh,norb_dz+1,lra-1)
      call link_b2_at_given_orb(mh)        !b^r
      call link_c1_to_given_orb(mh,lra+1,norb_inn)
      if (mh /= 0) then
        call value_sort_ploop(mh,.false.,.true.,.true.)
        if (log_prod == 3) then
          linelp = 16
          mhlp = mh
          nlg1 = lra
          nlg2 = 1
          call ext_head_in_dbl()
        else
          call save_lp(16,mh,lra,1)
          lpblock = lpblock+1
        end if
      end if
    end if

    logic_br(1) = .true.
    call link_c2_to_given_orb(mh,norb_dz+1,lra-1)
    call link_b2_at_given_orb(mh)        !b^r
    call link_c1_to_given_orb(mh,lra+1,norb_inn)
    if (mh /= 0) then
      call value_sort_ploop(mh,.false.,.true.,.true.)
      if (log_prod == 3) then
        linelp = 16
        mhlp = mh
        nlg1 = lra
        nlg2 = 2
        call ext_head_in_dbl()
      else
        call save_lp(16,mh,lra,2)
        lpblock = lpblock+1
      end if
    end if

  end if

  ! line=19 -d&r^l-
  call link_c1_to_given_orb(mh,norb_dz+1,lra-1)
  call link_d10_at_given_orb(mh)
  call link_c1_to_given_orb(mh,lra+1,norb_inn)
  if (mh == 0) cycle
  call value_sort_ploop(mh,.true.,.true.,.true.)
  if (log_prod == 3) then
    linelp = 19
    mhlp = mh
    nlg1 = lra
    nlg2 = 0
    call ext_head_in_dbl()
  else
    call save_lp(19,mh,lra,0)
    lpblock = lpblock+1
  end if
end do

do lrai=norb_dz+1,norb_inn-1
  lsmi = lsm_inn(lrai)
  do lraj=lrai+1,norb_inn
    lsmj = lsm_inn(lraj)
    lsmij = Mul(lsmi,lsmj)
    if (lsmij /= lsmact) cycle
    jk = ngw2(lrai-norb_frz)+ngw3(lraj-norb_frz)

    ! line=21 -b&l-b^r-
    call link_c1_to_given_orb(mh,norb_dz+1,lrai-1)
    call link_b3_at_given_orb(mh)        !b&l
    logic_br(1:mh) = .true.
    call link_c2_to_given_orb(mh,lrai+1,lraj-1)
    call link_b2_at_given_orb(mh)        !b^r
    call link_c1_to_given_orb(mh,lraj+1,norb_inn)
    if (mh == 0) cycle
    call value_sort_ploop(mh,.false.,.true.,.true.)
    if (log_prod == 3) then
      linelp = 21
      mhlp = mh
      nlg1 = jk
      nlg2 = 0
      call ext_head_in_dbl()
    else
      call save_lp(21,mh,jk,0)
      lpblock = lpblock+1
    end if
  end do
end do

return

end subroutine lp_head_in_dbl_2

!subroutine lp_head_in_dbl_2_mrpt2()          !for vd,dt,ds
!
!use gugaci_global, only: iml, imr, jml, jmr, jpadlr, jpadlrel, linelp, logic_br, lsm_inn, mhlp, ngw2, ngw3, nlg1, nlg2, norb_dz, &
!                         norb_frz, norb_inn
!use Symmetry_Info, only: Mul
!use Definitions, only: iwp
!
!implicit none
!integer(kind=iwp) :: imlr, jk, jmlr, lend, lra, lrai, lraj, lsma, lsmact, lsmi, lsmij, lsmj, lsta, mh
!logical(kind=iwp) :: do_16
!
!imlr = Mul(iml,imr)
!jmlr = Mul(jml,jmr)
!lsmact = Mul(imlr,jmlr)
!lsta = norb_dz+1
!lend = norb_inn
!
!do_16 = .false.
!select case (jpadlrel(jpadlr)) !200709
!  case default ! (1)
!    do lra=norb_dz+1,norb_inn
!      lsma = lsm_inn(lra)
!      ! line=19 -d&l^r-
!      call link_c1_to_given_orb(mh,norb_dz+1,lra-1)
!      call link_d10_at_given_orb(mh)
!      call link_c1_to_given_orb(mh,lra+1,norb_inn)
!      if (mh == 0) cycle
!      call value_sort_ploop(mh,.true.,.true.,.true.)
!      linelp = 19
!      mhlp = mh
!      nlg1 = lra
!      nlg2 = 0
!      !call print_lp()
!      call ext_head_in_dbl()
!    end do
!
!    do lrai=norb_dz+1,norb_inn-1
!      lsmi = lsm_inn(lrai)
!      do lraj=lrai+1,norb_inn
!        lsmj = lsm_inn(lraj)
!        lsmij = Mul(lsmi,lsmj)
!        if (lsmij /= lsmact) cycle
!        jk = ngw2(lrai-norb_frz)+ngw3(lraj-norb_frz)
!
!        ! line=21 -b&l-b^r-
!        call link_c1_to_given_orb(mh,norb_dz+1,lrai-1)
!        call link_b3_at_given_orb(mh)        !b&l
!        logic_br(1:mh) = .true.
!        call link_c2_to_given_orb(mh,lrai+1,lraj-1)
!        call link_b2_at_given_orb(mh)        !b^r
!        call link_c1_to_given_orb(mh,lraj+1,norb_inn)
!        if (mh == 0) cycle
!        call value_sort_ploop(mh,.false.,.true.,.true.)
!        linelp = 21
!        mhlp = mh
!        nlg1 = jk
!        nlg2 = 0
!        !call print_lp()
!        call ext_head_in_dbl()
!      end do
!    end do
!
!  case (2)
!    if (imlr == jmlr) then               !!! imlr == 1
!
!      ! line=13 -c'-
!      call link_c1_to_given_orb_coe(mh,lsta,lend)
!      if (mh == 0) then
!        do_16 = .true.
!      else
!        call value_sort_ploop(mh,.true.,.true.,.true.)
!        linelp = 13
!        mhlp = mh
!        nlg1 = 0
!        nlg2 = 0
!        !call print_lp()
!        call ext_head_in_dbl()
!      end if
!
!    end if
!
!  case (3)
!    do_16 = .true.
!
!  case (4)
!end select
!
!! line=16 -b^r-
!if (do_16) then
!  do lra=norb_dz+1,norb_inn
!    lsma = lsm_inn(lra)
!
!    if (lsma /= lsmact) cycle
!
!    if (jml == jmr) then
!      logic_br(1) = .false.
!      call link_c2_to_given_orb(mh,norb_dz+1,lra-1)
!      call link_b2_at_given_orb(mh)        !b^r
!      call link_c1_to_given_orb(mh,lra+1,norb_inn)
!      if (mh /= 0) then
!        call value_sort_ploop(mh,.false.,.true.,.true.)
!        linelp = 16
!        mhlp = mh
!        nlg1 = lra
!        nlg2 = 1
!        !call print_lp()
!        call ext_head_in_dbl()
!      end if
!    end if
!
!    logic_br(1) = .true.
!    call link_c2_to_given_orb(mh,norb_dz+1,lra-1)
!    call link_b2_at_given_orb(mh)        !b^r
!    call link_c1_to_given_orb(mh,lra+1,norb_inn)
!    if (mh == 0) cycle
!    call value_sort_ploop(mh,.false.,.true.,.true.)
!    linelp = 16
!    mhlp = mh
!    nlg1 = lra
!    nlg2 = 2
!    !call print_lp()
!    call ext_head_in_dbl()
!  end do
!end if
!
!return
!
!end subroutine lp_head_in_dbl_2_mrpt2

subroutine lp_head_in_dbl_3(ide)      !for ide=0:dd,tt,ide=1:ss,id

use gugaci_global, only: iml, imr, jml, jmr, logic_br, lpblock, lsm_inn, norb_dz, norb_inn
use Symmetry_Info, only: Mul
use Definitions, only: iwp

implicit none
integer(kind=iwp), intent(in) :: ide
integer(kind=iwp) :: imlr, jmlr, lend, lra, lsma, lsmact, lsta, mh

imlr = Mul(iml,imr)
jmlr = Mul(jml,jmr)
lsmact = Mul(imlr,jmlr)
lsta = norb_dz+1
lend = norb_inn

if (imlr == jmlr) then       !!! imlr == 1

  ! line=14 -c"-
  if (jml == jmr) then
    logic_br(1) = .false.
    call link_c2_to_given_orb(mh,lsta,lend)
    if (mh /= 0) then
      call value_sort_ploop(mh,.true.,.true.,.true.)
      call save_lp(14,mh,0,1)
      lpblock = lpblock+1
    end if
  end if

  logic_br(1) = .true.
  call link_c2_to_given_orb(mh,lsta,lend)
  if (mh /= 0) then
    call value_sort_ploop(mh,.true.,.true.,.true.)
    call save_lp(14,mh,0,2)
    lpblock = lpblock+1
  end if

end if

do lra=norb_dz+1,norb_inn
  lsma = lsm_inn(lra)
  if (lsma /= lsmact) cycle

  ! line=17 -b&l-
  call link_c1_to_given_orb(mh,norb_dz+1,lra-1)
  call link_b3_at_given_orb(mh)        !b&l
  logic_br(1:mh) = .true.
  call link_c2_to_given_orb(mh,lra+1,norb_inn)
  if (mh == 0) cycle
  select case (ide)
    case (0)
      call value_sort_ploop(mh,.false.,.true.,.true.)
    case default ! (1)
      call value_sort_ploop(mh,.false.,.true.,.false.)    !ss
    case (2)
      call value_sort_ploop(mh,.false.,.false.,.true.)    !st,ts
  end select
  call save_lp(17,mh,lra,0)
  lpblock = lpblock+1

end do

return

end subroutine lp_head_in_dbl_3

subroutine lp_head_in_dbl_4()

use gugaci_global, only: iml, imr, jml, jmr, logic_br, lpblock, lsm_inn, norb_dz, norb_inn
use Symmetry_Info, only: Mul
use Definitions, only: iwp

implicit none
integer(kind=iwp) :: imlr, jmlr, lend, lra, lsma, lsmact, lsta, mh

imlr = Mul(iml,imr)
jmlr = Mul(jml,jmr)
lsmact = Mul(imlr,jmlr)
lsta = norb_dz+1
lend = norb_inn

if (imlr == jmlr) then       !!! imlr == 1

  ! line=14 -c"-
  if (jml == jmr) then
    logic_br(1) = .false.
    call link_c2_to_given_orb(mh,lsta,lend)
    if (mh /= 0) then
      call value_sort_ploop(mh,.true.,.true.,.true.)
      call save_lp(14,mh,0,1)
      lpblock = lpblock+1
    end if
  end if

  logic_br(1) = .true.
  call link_c2_to_given_orb(mh,lsta,lend)
  if (mh /= 0) then
    call value_sort_ploop(mh,.true.,.true.,.true.)
    call save_lp(14,mh,0,2)
    lpblock = lpblock+1
  end if

end if
do lra=norb_dz+1,norb_inn
  lsma = lsm_inn(lra)
  if (lsma /= lsmact) cycle

  ! line=18 -b&r-
  call link_c1_to_given_orb(mh,norb_dz+1,lra-1)
  call link_b4_at_given_orb(mh)        !b&r
  logic_br(1:mh) = .true.
  call link_c2_to_given_orb(mh,lra+1,norb_inn)
  if (mh == 0) cycle
  call value_sort_ploop(mh,.false.,.true.,.true.)
  call save_lp(18,mh,lra,0)
  lpblock = lpblock+1

end do

return

end subroutine lp_head_in_dbl_4

!subroutine lp_head_in_dbl_4_mrpt2()
!
!use gugaci_global, only: iml, imr, jml, jmr, jpadlr, jpadlrel, linelp, logic_br, lsm_inn, mhlp, nlg1, nlg2, norb_dz, norb_inn
!use Symmetry_Info, only: Mul
!use Definitions, only: iwp
!
!implicit none
!integer(kind=iwp) :: imlr, jmlr, lend, lra, lsma, lsmact, lsta, mh
!logical(kind=iwp) :: do_15
!
!imlr = Mul(iml,imr)
!jmlr = Mul(jml,jmr)
!lsmact = Mul(imlr,jmlr)
!lsta = norb_dz+1
!lend = norb_inn
!
!do_15 = .false.
!select case (jpadlrel(jpadlr)) !200709
!  case default ! (1)
!    do_15 = .true.
!
!  case (2)
!
!  case (3)
!
!  case (4)
!    if (imlr == jmlr) then               !!! imlr == 1
!
!      ! line=14 -c"-
!      if (jml == jmr) then
!        logic_br(1) = .false.
!        call link_c2_to_given_orb(mh,lsta,lend)
!        if (mh /= 0) then
!          call value_sort_ploop(mh,.true.,.true.,.true.)
!          linelp = 14
!          mhlp = mh
!          nlg1 = 0
!          nlg2 = 1
!          call ext_head_in_dbl()
!        end if
!      end if
!
!      logic_br(1) = .true.
!      call link_c2_to_given_orb(mh,lsta,lend)
!      if (mh == 0) then
!        do_15 = .true.
!      else
!        call value_sort_ploop(mh,.true.,.true.,.true.)
!        linelp = 14
!        mhlp = mh
!        nlg1 = 0
!        nlg2 = 2
!        call ext_head_in_dbl()
!      end if
!
!    end if
!end select
!
!if (do_15) then
!  do lra=norb_dz+1,norb_inn
!    lsma = lsm_inn(lra)
!    if (lsma /= lsmact) cycle
!
!    ! line=18 -b&r-
!    call link_c1_to_given_orb(mh,norb_dz+1,lra-1)
!    call link_b4_at_given_orb(mh)        !b&r
!    logic_br(1:mh) = .true.
!    call link_c2_to_given_orb(mh,lra+1,norb_inn)
!    if (mh == 0) cycle
!    call value_sort_ploop(mh,.false.,.true.,.true.)
!    linelp = 18
!    mhlp = mh
!    nlg1 = lra
!    nlg2 = 0
!    call ext_head_in_dbl()
!
!  end do
!end if
!
!return
!
!end subroutine lp_head_in_dbl_4_mrpt2

!subroutine print_lp()
!
!use gugaci_global, only: ihy, ihyl, iml, imr, jml, jmr, jpadlr, jphy, jphyl, linelp, lp_count, lpnew_coe, lpnew_head, lpnew_lwei, &
!                         lpnew_rwei, mhlp, mhsum, mtype, nlg1, nlg2, norb_dz, norb_inn, ns_sm, nstaval, nvalue, vplpnew_w0, &
!                         vplpnew_w1
!use Symmetry_Info, only: Mul
!use Definitions, only: iwp
!
!implicit none
!integer(kind=iwp) :: in_, jhyl, jhyr, jph, jpml, jpmr, lg1, lg2, line, m, mh
!
!line = linelp
!mh = mhlp
!lg1 = nlg1
!lg2 = nlg2
!
!mhsum = mhsum+mh
!lp_count(line) = lp_count(line)+1
!if (((line == 14) .or. (line == 15) .or. (line == 16)) .and. (lg2 == 1)) then
!  mhsum = mhsum-mh
!  lp_count(line) = lp_count(line)-1
!end if
!!=======================================================================
!jpml = Mul(jml,ns_sm)
!jpmr = Mul(jmr,ns_sm)
!write(20,'(2x,10i8)') line,iml,imr,jpml,jpmr,jpadlr,mtype,mh,lg1,lg2
!write(20,'(7f12.6)') vplpnew_w0(1:mtype)
!write(20,'(7f12.6)') vplpnew_w1(1:mtype)
!write(20,'(8i4)') nstaval(1:mtype)
!write(20,'(8i4)') nvalue(1:mtype)
!write(20,'(8i6)') lpnew_lwei(1:mh)
!write(20,'(8i6)') lpnew_rwei(1:mh)
!!=======================================================================
!if (line <= 12) then
!  ! upwei_record
!  do m=1,mh
!    jph = lpnew_head(m)
!    jhyl = jphyl(jph)
!    jhyr = jphy(jph)
!    in_ = ihyl(jhyl)
!    write(20,'(8(1x,i4))') jph,jhyl,jhyr,in_,ihyl(jhyl+1:jhyl+in_),ihy(jhyr+1:jhyr+in_)
!  end do
!end if
!if ((line /= 1) .and. (line /= 13)) return
!! coe_record
!do m=1,mh
!  write(20,*) lpnew_coe(norb_dz+1:norb_inn,m)
!end do
!
!end subroutine print_lp

subroutine save_lp(line,mh,lg1,lg2)

use gugaci_global, only: idisk_lp, ihy, ihyl, iml, imr, jml, jmr, jpadlr, jphy, jphyl, lp_count, lpnew_coe, lpnew_head, &
                         lpnew_lwei, lpnew_rwei, LuLoop, mhlpmax, mhsum, mtype, norb_dz, norb_inn, ns_sm, nstaval, nvalue, &
                         vplpnew_w0, vplpnew_w1
use Symmetry_Info, only: Mul
use Definitions, only: iwp

implicit none
integer(kind=iwp), intent(in) :: line, mh, lg1, lg2
integer(kind=iwp) :: idum(1), in_, info(10), jhyl, jhyr, jph, jpml, jpmr, lenw, m

if (mh > mhlpmax) mhlpmax = mh
mhsum = mhsum+mh
lp_count(line) = lp_count(line)+1
if (((line == 14) .or. (line == 15) .or. (line == 16)) .and. (lg2 == 1)) then
  mhsum = mhsum-mh
  lp_count(line) = lp_count(line)-1
end if
info = 0
!=======================================================================
!jpml = Mul(jml,ns_sm)
!jpmr = Mul(jmr,ns_sm)
!write(200,'(2x,10i8)') line,iml,imr,jpml,jpmr,jpadlr,mtype,mh,lg1,lg2
!write(200,*) vplpnew_w0(1:mtype)
!write(200,*) vplpnew_w1(1:mtype)
!write(200,*) nstaval(1:mtype)
!write(200,*) nvalue(1:mtype)
!write(200,*) lpnew_lwei(1:mh)
!write(200,*) lpnew_rwei(1:mh)
!=======================================================================
jpml = Mul(jml,ns_sm)
jpmr = Mul(jmr,ns_sm)
info(1) = line
info(2) = iml
info(3) = imr
info(4) = jpml
info(5) = jpmr
info(6) = jpadlr
info(7) = mtype
info(8) = mh
info(9) = lg1
info(10) = lg2
call idafile(luloop,1,info,10,idisk_lp)
!write(u6,*) 'in save_lp, linelp',line,idisk_lp,mh,mtype
call ddafile(luloop,1,vplpnew_w0,mtype,idisk_lp)
call ddafile(luloop,1,vplpnew_w1,mtype,idisk_lp)
call idafile(luloop,1,nstaval,mtype,idisk_lp)
call idafile(luloop,1,nvalue,mtype,idisk_lp)
call idafile(luloop,1,lpnew_lwei,mh,idisk_lp)
call idafile(luloop,1,lpnew_rwei,mh,idisk_lp)
if (line <= 12) then
  do m=1,mh
    jph = lpnew_head(m)
    jhyl = jphyl(jph)
    jhyr = jphy(jph)
    in_ = ihyl(jhyl)
    idum(1) = in_
    call idafile(luloop,1,idum,1,idisk_lp)
    call idafile(luloop,1,ihyl(jhyl+1:jhyl+in_),in_,idisk_lp)
    call idafile(luloop,1,ihy(jhyr+1:jhyr+in_),in_,idisk_lp)
    !write(20) in,ihyl(jhyl+1:jhyl+in),ihy(jhyr+1:jhyr+in)
  end do
end if
if ((line /= 1) .and. (line /= 13)) return
!write(u6,*) 'in read_lp, write coe',idisk_lp,norb_inn-norb_dz
! coe_record
lenw = norb_inn-norb_dz
do m=1,mh
  call idafile(luloop,1,lpnew_coe(norb_dz+1:norb_inn,m),lenw,idisk_lp)
end do

return

end subroutine save_lp

subroutine read_lp()

use gugaci_global, only: idisk_lp, ihy, ihyl, iml, imr, jml, jmr, jpadlr, jphy, linelp, lpnew_coe, lpnew_lwei, lpnew_rwei, LuLoop, &
                         mhlp, mtype, ndim, nlg1, nlg2, norb_dz, norb_inn, nstaval, nvalue, vplpnew_w0, vplpnew_w1
use Definitions, only: iwp

implicit none
integer(kind=iwp) :: ihypos, info(10), lenr, m

call idafile(luloop,2,info,10,idisk_lp)
linelp = info(1)
iml = info(2)
imr = info(3)
jml = info(4)
jmr = info(5)
jpadlr = info(6)
mtype = info(7)
mhlp = info(8)
nlg1 = info(9)
nlg2 = info(10)
call ddafile(luloop,2,vplpnew_w0,mtype,idisk_lp)
call ddafile(luloop,2,vplpnew_w1,mtype,idisk_lp)
call idafile(luloop,2,nstaval,mtype,idisk_lp)
call idafile(luloop,2,nvalue,mtype,idisk_lp)
call idafile(luloop,2,lpnew_lwei,mhlp,idisk_lp)
call idafile(luloop,2,lpnew_rwei,mhlp,idisk_lp)
!=======================================================================

if (linelp <= 12) then
  ! upwei_read
  ihypos = 1
  do m=1,mhlp
    jphy(m) = ihypos
    call idafile(luloop,2,info,1,idisk_lp)
    ndim = info(1)
    call idafile(luloop,2,ihyl(ihypos+1:ihypos+ndim),ndim,idisk_lp)
    call idafile(luloop,2,ihy(ihypos+1:ihypos+ndim),ndim,idisk_lp)
    ihy(ihypos) = ndim
    ihypos = ihypos+ndim+1
  end do
end if
if ((linelp /= 1) .and. (linelp /= 13)) return
!write(u6,*) 'in read_lp, read coe',idisk_lp,norb_inn-norb_dz
! coe_read
lenr = norb_inn-norb_dz
do m=1,mhlp
  call idafile(luloop,2,lpnew_coe(norb_dz+1:norb_inn,m),lenr,idisk_lp)
end do

return

end subroutine read_lp

subroutine value_sort_ploop(mh,logic_ar,logic_w0,logic_w1)

use gugaci_global, only: lp_coe, lp_head, lp_lwei, lp_rwei, lpnew_coe, lpnew_head, lpnew_lwei, lpnew_rwei, mtype, norb_dz, &
                         norb_inn, nstaval, nvaltype, nvalue, vplp_w0, vplp_w1, vplpnew_w0, vplpnew_w1
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: mh
logical(kind=iwp), intent(in) :: logic_ar, logic_w0, logic_w1
integer(kind=iwp) :: lr, m, mnew, mty, n, nty
real(kind=wp) :: s_w0, s_w1, w0, w1
integer(kind=iwp), allocatable :: ntype(:)

call mma_allocate(ntype,nvaltype,label='nvaltype')
nvalue(1:nvaltype) = 0
ntype(1:mh) = 0
mtype = 0
if (.not. logic_w0) then
  do n=1,mh
    vplp_w0(n) = Zero
  end do
end if
if (.not. logic_w1) then
  do n=1,mh
    vplp_w1(n) = Zero
  end do
end if
do n=1,mh
  if (ntype(n) /= 0) cycle
  s_w0 = vplp_w0(n)
  s_w1 = vplp_w1(n)
  mtype = mtype+1
  if (mtype > nvaltype) then
    write(u6,*) ' out of array boundary'
    write(u6,*) ' subroutine value_sort_ploop'
    write(u6,*) ' program stop'
#   ifndef MOLPRO
    call abend()
#   endif
  end if
  ntype(n) = mtype
  nvalue(mtype) = nvalue(mtype)+1
  vplpnew_w0(mtype) = vplp_w0(n)
  vplpnew_w1(mtype) = vplp_w1(n)
  do m=n+1,mh
    if (ntype(m) /= 0) cycle
    w0 = vplp_w0(m)
    w1 = vplp_w1(m)
    if ((w0 /= s_w0) .or. (w1 /= s_w1)) cycle
    ntype(m) = mtype
    nvalue(mtype) = nvalue(mtype)+1
  end do
end do

nstaval(1) = 0
do mty=1,mtype-1
  nstaval(mty+1) = nstaval(mty)+nvalue(mty)
end do

mnew = 0
do mty=1,mtype
  do n=1,mh
    nty = ntype(n)
    if (nty == mty) then
      mnew = mnew+1
      lpnew_head(mnew) = lp_head(n)
      lpnew_lwei(mnew) = lp_lwei(n)
      lpnew_rwei(mnew) = lp_rwei(n)
      if (logic_ar) then
        do lr=norb_dz+1,norb_inn
          lpnew_coe(lr,mnew) = lp_coe(lr,n)
        end do
      end if
    end if
  end do
end do
call mma_deallocate(ntype)

return

end subroutine value_sort_ploop

subroutine ext_head_in_act()

use gugaci_global, only: ipaety, jml, jmr, logic_dh, ns_sm
use Symmetry_Info, only: Mul

implicit none

logic_dh = .false.
jml = Mul(jml,ns_sm)
jmr = Mul(jmr,ns_sm)
select case (ipaety)
  case default ! (10)
    call sv_ext_head_in_act()
  case (17)
    call tv_ext_head_in_act()
  case (23)
    call dv_ext_head_in_act()
  case (26)
    call vd_ext_head_in_act()
  case (1:9,11:16,18:22,24:25)
end select
jml = Mul(jml,ns_sm)
jmr = Mul(jmr,ns_sm)

return

end subroutine ext_head_in_act

subroutine ext_head_in_dbl()

use gugaci_global, only: ipaety, jml, jmr, logic_dh, ns_sm
use Symmetry_Info, only: Mul

implicit none

logic_dh = .true.
jml = Mul(jml,ns_sm)
jmr = Mul(jmr,ns_sm)
select case (ipaety)
  case default ! (10)
    call sv_ext_head_in_dbl()
  case (17)
    call tv_ext_head_in_dbl()
  case (23)
    call dv_ext_head_in_dbl()
  case (26)
    call vd_ext_head_in_dbl()
  case (1:9,11:16,18:22,24:25)
end select
jml = Mul(jml,ns_sm)
jmr = Mul(jmr,ns_sm)

return

end subroutine ext_head_in_dbl
