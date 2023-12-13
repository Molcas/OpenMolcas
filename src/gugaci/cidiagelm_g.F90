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

subroutine diagonal_loop_wyb_g()  !  for norb_act<>0

use gugaci_global, only: ipae, iseg_downwei, iw_downwei, iw_sta, jd, jpad, jpad_upwei, jpae, js, jt, jv, mxnode, ndim, ng_sm, &
                         nu_ad, nu_ae
use Definitions, only: iwp, iwp

implicit none
integer(kind=iwp) :: im, iwupwei, jaedownwei, jpad_, ndimsum

!do lr0=2,norb_all
!  do lr=1,lr0-1
!    vdint(lr,lr0) = voint(lr0,lr)-vdint(lr0,lr)-vdint(lr0,lr)   ! 520
!    write(u6,'(2i4,3f14.8)') lr,lr0,voint(lr0,lr),vdint(lr0,lr),vdint(lr,lr0)
!  end do
!end do
!write(u6,*) '               ***** start h-diaelm *****'
!**************lyb***************
!     vector1(1:nci_dim)=vpotnuc
!********************************
!wl8 = hnil*(hnil-1)*vmd(lr,lr)*Half+hnil*vo(lr,lr)   ! 800
!write(u6,*) '         jpad,     jpae,     ndim,      nohy'

ndimsum = 0
jpae = jv
ipae = 1
jaedownwei = iseg_downwei(ipae)
do jpad_=1,mxnode
  jpad = jpad_ ! jpad is in global module, is this necessary?
  iw_sta(jpad,ipae) = ndimsum
  if (nu_ad(jpad) == 0) cycle
  call seg_drt()
  iwupwei = jpad_upwei(jpad)
  iw_downwei(jpad,ipae) = ndim
  ndimsum = ndimsum+ndim*jaedownwei*iwupwei
  if (ndim == 0) cycle
  call diagonal_act_d_g()
  call diagonal_act_c_g()
end do
do im=1,ng_sm
  jpae = jd(im)
  ipae = 1+im
  if (nu_ae(ipae) == 0) cycle
  jaedownwei = iseg_downwei(ipae)
  do jpad_=1,mxnode
    jpad = jpad_ ! jpad is in global module, is this necessary?
    iw_sta(jpad,ipae) = ndimsum
    if (nu_ad(jpad) == 0) cycle
    call seg_drt()
    iwupwei = jpad_upwei(jpad)
    iw_downwei(jpad,ipae) = ndim
    !if (jpad >= 26) write(u6,*)
    ndimsum = ndimsum+ndim*jaedownwei*iwupwei
    if (ndim == 0) cycle
    call diagonal_act_d_g()
    call diagonal_act_c_g()
  end do
end do

do im=1,ng_sm
  jpae = jt(im)
  ipae = 9+im
  if (nu_ae(ipae) == 0) cycle
  jaedownwei = iseg_downwei(ipae)
  do jpad_=1,mxnode
    jpad = jpad_ ! jpad is in global module, is this necessary?
    iw_sta(jpad,ipae) = ndimsum
    if (nu_ad(jpad) == 0) cycle
    call seg_drt()
    iwupwei = jpad_upwei(jpad)
    iw_downwei(jpad,ipae) = ndim
    ndimsum = ndimsum+ndim*jaedownwei*iwupwei
    if (ndim == 0) cycle
    call diagonal_act_d_g()
    call diagonal_act_c_g()
  end do
end do
do im=1,ng_sm
  jpae = js(im)
  ipae = 17+im
  jaedownwei = iseg_downwei(ipae)
  if (nu_ae(ipae) == 0) cycle
  do jpad_=1,mxnode
    jpad = jpad_ ! jpad is in global module, is this necessary?
    iw_sta(jpad,ipae) = ndimsum
    if (nu_ad(jpad) == 0) cycle
    call seg_drt()
    iwupwei = jpad_upwei(jpad)
    iw_downwei(jpad,ipae) = ndim
    ndimsum = ndimsum+ndim*jaedownwei*iwupwei
    if (ndim == 0) cycle
    call diagonal_act_d_g()
    call diagonal_act_c_g()
  end do
end do
call diagonal_dbl_g()
call diagonal_ext_g()

return

end subroutine diagonal_loop_wyb_g

subroutine diagonal_act_c_g()

use gugaci_global, only: iy, jb, jeh, jj_sub, jpad, jph, jwh, maxpl, norb_act, norb_dz, norb_inn, th, thh
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp) :: idl, ind1, isq, iwa, je, jeb, jp, jpb, lr, m, me, mh, mp, mpe, mw
real(kind=wp) :: vlop0, vlop1, w, ww
integer(kind=iwp), allocatable :: jee(:), jpe(:), jwe(:)
real(kind=wp), allocatable :: te(:), tee(:)

!write(u6,*) '               ***** start h-diaelm *****'
!write(u6,*) jpad,jpae
!ndr(:) = 0
if (norb_act == 0) then
  mh = 1
  th(1) = One
  thh(1) = One
  call diagonal_link_dae_g(mh)
  return
end if
mp = 0
! 520
mh = 0
me = 0
jp = jpad
jpb = jb(jp)
do idl=1,4
  if (jj_sub(idl,jp) == 0) cycle
  ind1 = idl
  ! link c"
  call smidc2(isq,w,ww,mw,ind1,jpb)
  mh = mh+1
  jeh(mh) = jj_sub(idl,jp)
  th(mh) = w
  thh(mh) = ww
  jph(mh) = 0
  jwh(mh) = 0
  if (idl /= 1) jwh(mh) = iy(idl,jp)
  ! complete a loop 'v'
  if (ind1 == 1) cycle
  call stml(isq,w,ww,mw,ind1-1,jpb)
  vlop0 = w
  vlop1 = ww
  if ((vlop0 == Zero) .and. (vlop1 == Zero)) cycle
  mpe = jj_sub(idl,jp)
  iwa = iy(idl,jp)
  call diagonal_link_ad_g(mpe,iwa,vlop0,vlop1)
end do
!***********************************************************************
! write(u6,*) ad(i)
!***********************************************************************
call mma_allocate(te,maxpl,label='te')
call mma_allocate(tee,maxpl,label='tee')
call mma_allocate(jpe,maxpl,label='jpe')
call mma_allocate(jee,maxpl,label='jee')
call mma_allocate(jwe,maxpl,label='jwe')
lr = norb_dz+1
do
  if (lr == norb_inn) then
    if (mh /= 0) call diagonal_link_dae_g(mh)
    exit
  end if
  lr = lr+1
  me = 0
  do m=1,mh
    je = jeh(m)
    jeb = jb(je)
    !jp = jph(m)
    do idl=1,4
      if (jj_sub(idl,je) == 0) cycle
      ind1 = idl
      !if (lr /= 1) then
      ! link loop
      call smidc2(isq,w,ww,mw,ind1,jeb)
      me = me+1
      jwe(me) = jwh(m)
      if (idl /= 1) jwe(me) = jwe(me)+iy(idl,je)
      jee(me) = jj_sub(idl,je)
      te(me) = th(m)*w
      tee(me) = thh(m)*ww
      jpe(me) = jph(m)
      !end if
      ! complete a loop 'v'
      if (ind1 == 1) cycle
      call stml(isq,w,ww,mw,ind1-1,jeb)
      vlop0 = th(m)*w
      vlop1 = thh(m)*ww
      if ((vlop0 == Zero) .and. (vlop1 == Zero)) cycle
      mp = mp+1
      mpe = jj_sub(idl,je)
      iwa = jwh(m)
      if (idl /= 1) iwa = iwa+iy(idl,je)
      call diagonal_link_ad_g(mpe,iwa,vlop0,vlop1)
      !*****   520  ****************************************************
    end do
  end do
  !*********************************************************************
  ! write(u6,*) ad(i)
  !*********************************************************************
  do m=1,me
    th(m) = te(m)
    te(m) = Zero
    thh(m) = tee(m)
    tee(m) = Zero
    jwh(m) = jwe(m)
    jwe(m) = 0
    jeh(m) = jee(m)
    jee(m) = 0
    jph(m) = jpe(m)
    jpe(m) = 0
  end do
  mh = me
  !if (ndr(lr) < mh) ndr(lr) = mh
end do
!do m=1,mh
!  th(m) = Zero
!  thh(m) = Zero
!  jwh(m) = 0
!  jph(m) = 0
!  jeh(m) = 0
!end do

call mma_deallocate(te)
call mma_deallocate(tee)
call mma_deallocate(jpe)
call mma_deallocate(jee)
call mma_deallocate(jwe)

return

end subroutine diagonal_act_c_g

subroutine diagonal_act_d_g()

use gugaci_global, only: iy, jb, jeh, jj_sub, jph, jwh, maxpl, no, norb_dz, norb_inn, th, thh
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One, Two
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp) :: idl, ind1, isq, iwa, jbl, je, jeb, jp, jp0, jp1, jw, lr, lr0, m, me, mh, mp, mpe, mw, nxo
real(kind=wp) :: vlop0, vlop1, w, wt, ww
integer(kind=iwp), allocatable :: jee(:), jpe(:), jwe(:)
real(kind=wp), allocatable :: te(:), tee(:)

!write(u6,*) '               ***** start h-diaelm *****'
!write(u6,*) '   diagonal_act_d:',jpad,ipae
!ndr(:) = 0
do lr=norb_dz+1,norb_inn
  jp0 = no(lr-1)+1
  jp1 = no(lr)
  lr0 = lr
  do jp=jp0,jp1
    if (iy(1,jp) == 0) cycle
    ! 800
    do idl=2,3
      mpe = jj_sub(idl,jp)
      if (mpe == 0) cycle
      !wt = voint(lr0,lr0)    ! hnil=1
      wt = One
      jw = iy(idl,jp)
      call prodel_1(3,wt,jp,mpe,jw,lr0,lr0)
    end do
    mpe = jj_sub(4,jp)
    if (mpe /= 0) then
      !wt = vdint(lr0,lr0)+Two*voint(lr0,lr0)     !idl=4 hnil=2
      jw = iy(4,jp)
      wt = Two
      call prodel_1(3,wt,jp,mpe,jw,lr0,lr0)
      wt = One
      call trans_ijkl_intpos(lr0,lr0,lr0,lr0,nxo)
      call prodel_2(3,wt,jp,mpe,jw,nxo)
    end if
  end do
end do
call mma_allocate(te,maxpl,label='te')
call mma_allocate(tee,maxpl,label='tee')
call mma_allocate(jpe,maxpl,label='jpe')
call mma_allocate(jee,maxpl,label='jee')
call mma_allocate(jwe,maxpl,label='jwe')
! 520
do lr0=norb_dz+1,norb_inn
  mh = 0
  me = 0
  jp0 = no(lr0-1)+1
  jp1 = no(lr0)
  if (jp0 > jp1) jp0 = jp1
  do jp=jp0,jp1
    if (iy(1,jp) == 0) cycle
    !start '^'
    do idl=1,4
      if (jj_sub(idl,jp) == 0) cycle
      jbl = jb(jp)
      ind1 = idl-1
      if (ind1 == 0) cycle
      call stmh(isq,w,ww,mw,ind1,jbl)
      mh = mh+1
      jeh(mh) = jj_sub(idl,jp)
      th(mh) = w
      thh(mh) = ww
      jph(mh) = jp
      jwh(mh) = 0
      if (idl /= 1) jwh(mh) = iy(idl,jp)
    end do
  end do
  !*********************************************************************
  !write(u6,*) ad(i)
  !*********************************************************************
  lr = lr0
  !if (ndr(lr) < mh) ndr(lr) = mh
  do
    if (lr == norb_inn) then
      call diagonal_link_ae_g(mh)
      exit
    end if
    lr = lr+1
    me = 0

    mp = 0
    do m=1,mh
      je = jeh(m)
      jeb = jb(je)
      jp = jph(m)
      do idl=1,4
        if (jj_sub(idl,je) == 0) cycle
        ind1 = idl
        if (lr /= 1) then
          ! link loop
          call smidc2(isq,w,ww,mw,ind1,jeb)
          me = me+1
          jwe(me) = jwh(m)
          if (idl /= 1) jwe(me) = jwe(me)+iy(idl,je)
          jee(me) = jj_sub(idl,je)
          te(me) = th(m)*w
          tee(me) = thh(m)*ww
          jpe(me) = jph(m)
        end if
        ! complete a loop 'v'
        if (ind1 == 1) cycle
        call stml(isq,w,ww,mw,ind1-1,jeb)
        vlop0 = th(m)*w
        vlop1 = thh(m)*ww
        if ((vlop0 == Zero) .and. (vlop1 == Zero)) cycle
        mp = mp+1
        mpe = jj_sub(idl,je)
        iwa = jwh(m)
        if (idl /= 1) iwa = iy(idl,je)+iwa
        !wt = (vlop0-vlop1)*voint(lr,lr0)-Two*vlop0*vdint(lr,lr0)
        wt = vlop0-vlop1
        call trans_ijkl_intpos(lr,lr0,lr,lr0,nxo)
        !write(nf2,'(4i4,i8)') lr,lr0,lr,lr0,nxo
        call prodel_2(3,wt,jp,mpe,iwa,nxo)
        wt = -Two*vlop0
        call trans_ijkl_intpos(lr,lr,lr0,lr0,nxo)
        call prodel_2(3,wt,jp,mpe,iwa,nxo)

        !*****   520  **************************************************
      end do
    end do
    !*******************************************************************
    !write(u6,*) ad(i)
    !*******************************************************************
    do m=1,me
      th(m) = te(m)
      te(m) = Zero
      thh(m) = tee(m)
      tee(m) = Zero
      jwh(m) = jwe(m)
      jwe(m) = 0
      jeh(m) = jee(m)
      jee(m) = 0
      jph(m) = jpe(m)
      jpe(m) = 0
    end do
    mh = me
    !if (ndr(lr) < mh) ndr(lr) = mh
  end do
  !do m=1,mh
  !  th(m) = Zero
  !  thh(m) = Zero
  !  jwh(m) = 0
  !  jph(m) = 0
  !  jeh(m) = 0
  !end do
end do

call mma_deallocate(te)
call mma_deallocate(tee)
call mma_deallocate(jpe)
call mma_deallocate(jee)
call mma_deallocate(jwe)

return

end subroutine diagonal_act_d_g

subroutine diagonal_link_ae_g(mh)

use gugaci_global, only: ibsm_ext, iesm_ext, ipae, jph, jwh, kk, lsm, ng_sm, nlsm_ext, norb_all, norb_ext, th, thh, v_onevsqtwo, &
                         v_sqthreevsqtwo, v_sqtwo
use Symmetry_Info, only: Mul
use Constants, only: Two
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(inout) :: mh
integer(kind=iwp) :: ima, imae, imb, ip, ityae, iwa, iwe, jp, la, lb, lbend, lbsta, lr0, lra, lrb, ma, nxo
real(kind=wp) :: vlop0, vlop1, wg13, wg14, wg38, wg50, wld, wls, wlt, wwg38, wwg50

if (ipae == 1) return   !could not link
ityae = (ipae-1)/8+1
imae = mod(ipae-1,8)
if (imae == 0) then
  ityae = ityae-1
  imae = 8
end if
do ip=1,mh
  iwe = 0
  jp = jph(ip)
  lr0 = kk(jp)
  !write(u6,*) 'ip,jpe,ind0',ip,jpe,ind0
  iwa = jwh(ip)
  vlop0 = th(ip)
  vlop1 = thh(ip)

  !wl5 = (vlop0-vlop1)*vo(lr0,lr)-Two*vlop0*vmd(lr0,lr)
  !wl8 = vlop0*(vo(lr0,lr0)+(vlop0-1)*Half*vmd(lr0,lr0))
  ! two-index,one-loop 520
  ! 520 = <a,j,k,a>:13,14(ss=3),38(tt=2),50(dd=1)
  !write(nf2,*) 'ityae',ityae

  !link arc_d
  select case (ityae)
    case default ! (1)
      !zz = '  g50  '
      wg50 = vlop0*v_onevsqtwo
      wwg50 = -vlop1*v_sqthreevsqtwo
      do la=1,norb_ext
        ma = lsm(la)
        if (ma /= imae) cycle
        lra = norb_all-la+1
        iwe = iwe+1
        !wld = (wg50-wwg50)*voint(lra,lr0)-Two*wg50*vdint(lra,lr0)
        wld = wg50-wwg50
        call trans_ijkl_intpos(lra,lr0,lra,lr0,nxo)
        call prodel_2(5,wld,jp,iwa,iwe,nxo)
        wld = -Two*wg50
        call trans_ijkl_intpos(lra,lra,lr0,lr0,nxo)
        call prodel_2(5,wld,jp,iwa,iwe,nxo)

        !write(u6,'(a11,2i3,i6,1x,5f10.4)') zz,lr0,la,jwl,vo(lr0,la),vmd(lr0,la),wg50,wwg50,wl
      end do

    case (2)
      !zz = '  g38,39  '
      wg38 = -vlop0*v_onevsqtwo
      wwg38 = vlop1
      do ima=1,ng_sm
        imb = Mul(ima,imae)
        if (imb > ima) cycle
        do la=ibsm_ext(ima),iesm_ext(ima)
          lra = norb_all-la+1
          lbsta = ibsm_ext(imb)
          lbend = iesm_ext(imb)
          if (ima == imb) lbend = la-1
          do lb=lbsta,lbend
            lrb = norb_all-lb+1
            iwe = iwe+1
            !wlt = (wg38-wwg38)*(voint(lra,lr0)+voint(lrb,lr0))-Two*wg38*(vdint(lra,lr0)+vdint(lrb,lr0))
            wlt = wg38-wwg38
            call trans_ijkl_intpos(lra,lr0,lra,lr0,nxo)
            call prodel_2(5,wlt,jp,iwa,iwe,nxo)
            wlt = wg38-wwg38
            call trans_ijkl_intpos(lrb,lr0,lrb,lr0,nxo)
            call prodel_2(5,wlt,jp,iwa,iwe,nxo)
            wlt = -Two*wg38
            call trans_ijkl_intpos(lra,lra,lr0,lr0,nxo)
            call prodel_2(5,wlt,jp,iwa,iwe,nxo)
            wlt = -Two*wg38
            call trans_ijkl_intpos(lrb,lrb,lr0,lr0,nxo)
            call prodel_2(5,wlt,jp,iwa,iwe,nxo)
          end do
        end do
      end do

    case (3)
      !zz = '  g14,15  '
      wg14 = -vlop0*v_onevsqtwo
      do ima=1,ng_sm
        imb = Mul(ima,imae)
        if (imb > ima) cycle
        if (nlsm_ext(ima) == 0) cycle
        if (nlsm_ext(imb) == 0) cycle
        do la=ibsm_ext(ima),iesm_ext(ima)
          lra = norb_all-la+1
          lbsta = ibsm_ext(imb)
          lbend = iesm_ext(imb)
          if (ima == imb) lbend = la-1
          do lb=lbsta,lbend
            lrb = norb_all-lb+1
            iwe = iwe+1
            !wls = wg14*(voint(lra,lr0)+voint(lrb,lr0))-Two*wg14*(vdint(lra,lr0)+vdint(lrb,lr0))
            wls = wg14
            call trans_ijkl_intpos(lra,lr0,lra,lr0,nxo)
            call prodel_2(5,wls,jp,iwa,iwe,nxo)
            wls = wg14
            call trans_ijkl_intpos(lrb,lr0,lrb,lr0,nxo)
            call prodel_2(5,wls,jp,iwa,iwe,nxo)
            wls = -Two*wg14
            call trans_ijkl_intpos(lra,lra,lr0,lr0,nxo)
            call prodel_2(5,wls,jp,iwa,iwe,nxo)
            wls = -Two*wg14
            call trans_ijkl_intpos(lrb,lrb,lr0,lr0,nxo)
            call prodel_2(5,wls,jp,iwa,iwe,nxo)
          end do
        end do
      end do
      if (ipae /= 18) cycle
      !zz = '  g13     '
      wg13 = -vlop0*v_sqtwo
      do la=1,norb_ext
        lra = norb_all-la+1
        iwe = iwe+1
        !wls = wg13*(voint(lra,lr0)-Two*vdint(lra,lr0))
        wls = wg13
        call trans_ijkl_intpos(lra,lr0,lra,lr0,nxo)
        call prodel_2(5,wls,jp,iwa,iwe,nxo)
        wls = -Two*wg13
        call trans_ijkl_intpos(lra,lra,lr0,lr0,nxo)
        call prodel_2(5,wls,jp,iwa,iwe,nxo)

        !write(u6,*) ' g13 ',vo(lr0,la),vo(lr0,lb),vmd(lr0,la),vmd(lr0,lb)
      end do
  end select
end do
mh = 0

return

end subroutine diagonal_link_ae_g

subroutine diagonal_link_ad_g(mpe,iwa,vlop0,vlop1)

use gugaci_global, only: fg, jb_sys, jpad, jud, just, kk, lsm_inn, norb_dz, norb_frz, ns_sm, pd, pdd, ps1, ps2, ps3, ps4, pt, ptt, &
                         v_onevsqtwo, v_sqtwo
use Symmetry_Info, only: Mul
use Constants, only: Zero, Two
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: mpe, iwa
real(kind=wp), intent(in) :: vlop0, vlop1
integer(kind=iwp) :: imad, imd, imi, imij, imj, ityad, iwd, iws, iwt, lr, lra, lri, lrj, lrjsta, nxo
real(kind=wp) :: vl0, vl1, wld, wls, wlt, wlv
real(kind=wp) :: fqi

if (norb_dz == 0) return
ityad = 1
if (jpad /= 1) ityad = (jpad-1)/8+2
imad = mod(jpad-1,8)
if (imad == 0) then
  ityad = ityad-1
  imad = 8
end if
lra = kk(mpe)-1
select case (ityad)
  case default ! (1)
    ! v: d&r&l(3)
    fqi = -fg
    vl0 = fqi*v_sqtwo*vlop0
    wlv = Zero
    do lr=1,norb_dz
      !wlv = wlv-vl0*vdint(lr,lra)
      wlv = -vl0
      call trans_ijkl_intpos(lra,lr,lra,lr,nxo)
      !write(nf2,'(4i4,i8)') lra,lr,lra,lr,nxo
      call prodel_2(4,wlv,mpe,0,iwa,nxo)
      wlv = Two*vl0
      call trans_ijkl_intpos(lra,lra,lr,lr,nxo)
      call prodel_2(4,wlv,mpe,0,iwa,nxo)
    end do
    !call prodel(4,wlv,mpe,0,iwa)

  case (2)
    !jpad = jd(im)
    fqi = -fg
    do lri=norb_frz+1,norb_dz
      imd = Mul(lsm_inn(lri),ns_sm)
      if (imd /= imad) cycle
      iwd = jud(lri)

      ! d: d&r&l(2)
      vl0 = fqi*v_onevsqtwo*vlop0
      vl1 = pd*vlop1
      !wld = -Two*vl0*vdint(lra,lri)+(vl0-vl1)*voint(lra,lri)
      wld = -Two*vl0
      call trans_ijkl_intpos(lra,lra,lri,lri,nxo)
      call prodel_2(4,wld,mpe,iwd,iwa,nxo)
      wld = vl0-vl1
      call trans_ijkl_intpos(lra,lri,lra,lri,nxo)
      call prodel_2(4,wld,mpe,iwd,iwa,nxo)

      ! d: d&r&l(3)+c"(2)
      vl0 = fqi*v_sqtwo*vlop0
      do lr=1,lri-1
        !wld = wld+vl0*(voint(lra,lr)-Two*vdint(lra,lr))
        wld = -Two*vl0
        call trans_ijkl_intpos(lra,lra,lr,lr,nxo)
        call prodel_2(4,wld,mpe,iwd,iwa,nxo)
        wld = vl0
        call trans_ijkl_intpos(lra,lr,lra,lr,nxo)
        call prodel_2(4,wld,mpe,iwd,iwa,nxo)
      end do
      ! d: d&r&l(3)
      do lr=lri+1,norb_dz
        !wld = wld+vl0*(voint(lra,lr)-Two*vdint(lra,lr))
        wld = -Two*vl0
        call trans_ijkl_intpos(lra,lra,lr,lr,nxo)
        call prodel_2(4,wld,mpe,iwd,iwa,nxo)
        wld = vl0
        call trans_ijkl_intpos(lra,lr,lra,lr,nxo)
        call prodel_2(4,wld,mpe,iwd,iwa,nxo)
      end do

      !call prodel(4,wld,mpe,iwd,iwa)
    end do

  case (3)
    !jpad = jt(im)
    fqi = fg
    iwt = 0
    do lri=norb_frz+1,norb_dz
      imi = Mul(lsm_inn(lri),ns_sm)
      do lrj=lri+1,norb_dz
        imj = lsm_inn(lrj)
        imij = Mul(imi,imj)
        if (imij /= imad) cycle
        iwt = just(lri,lrj)

        ! t: d&r&l(2)
        ! t: d&r&l(2)+c"(2)
        vl0 = fqi*v_onevsqtwo*vlop0
        vl1 = pt*vlop1
        !wlt = -Two*vl0*vdint(lra,lri)+(vl0-vl1)*voint(lra,lri)
        wlt = -Two*vl0
        call trans_ijkl_intpos(lra,lra,lri,lri,nxo)
        call prodel_2(4,wlt,mpe,iwt,iwa,nxo)
        wlt = vl0-vl1
        call trans_ijkl_intpos(lra,lri,lra,lri,nxo)
        call prodel_2(4,wlt,mpe,iwt,iwa,nxo)

        !wlt = wlt-Two*vl0*vdint(lra,lrj)+(vl0-vl1)*voint(lra,lrj)
        wlt = -Two*vl0
        call trans_ijkl_intpos(lra,lra,lrj,lrj,nxo)
        call prodel_2(4,wlt,mpe,iwt,iwa,nxo)
        wlt = vl0-vl1
        call trans_ijkl_intpos(lra,lrj,lra,lrj,nxo)
        call prodel_2(4,wlt,mpe,iwt,iwa,nxo)

        vl0 = fqi*v_sqtwo*vlop0
        ! t: d&r&l(3)+c"(2)+c"(2)
        do lr=1,lri-1
          !wlt = wlt+vl0*vdint(lr,lra)
          wlt = vl0
          call trans_ijkl_intpos(lra,lr,lra,lr,nxo)
          call prodel_2(4,wlt,mpe,iwt,iwa,nxo)
          wlt = -Two*vl0
          call trans_ijkl_intpos(lra,lra,lr,lr,nxo)
          call prodel_2(4,wlt,mpe,iwt,iwa,nxo)
        end do
        ! t: d&r&l(3)+c"(2)
        do lr=lri+1,lrj-1
          !wlt = wlt+vl0*vdint(lr,lra)
          wlt = vl0
          call trans_ijkl_intpos(lra,lr,lra,lr,nxo)
          call prodel_2(4,wlt,mpe,iwt,iwa,nxo)
          wlt = -Two*vl0
          call trans_ijkl_intpos(lra,lra,lr,lr,nxo)
          call prodel_2(4,wlt,mpe,iwt,iwa,nxo)
        end do
        ! t: d&r&l(3)
        do lr=lrj+1,norb_dz
          !wlt = wlt+vl0*vdint(lr,lra)
          wlt = vl0
          call trans_ijkl_intpos(lra,lr,lra,lr,nxo)
          call prodel_2(4,wlt,mpe,iwt,iwa,nxo)
          wlt = -Two*vl0
          call trans_ijkl_intpos(lra,lra,lr,lr,nxo)
          call prodel_2(4,wlt,mpe,iwt,iwa,nxo)
        end do
        !call prodel(4,wlt,mpe,iwt,iwa)
      end do
    end do

  case (4)
    fqi = fg
    iws = 0
    do lri=norb_frz+1,norb_dz
      if (imad == ns_sm) then
        lrj = lri
        iws = just(lri,lri)
        ! s: d&r&l(3)+c"(0)
        ! s: d&r&l(3)
        vl0 = fqi*v_sqtwo*vlop0
        !wls = Zero
        do lr=1,norb_dz
          if (lr == lri) cycle
          !wls = wls+vl0*(voint(lra,lr)-Two*vdint(lra,lr))
          wls = vl0
          call trans_ijkl_intpos(lra,lr,lra,lr,nxo)
          call prodel_2(4,wls,mpe,iws,iwa,nxo)
          wls = -Two*vl0
          call trans_ijkl_intpos(lra,lra,lr,lr,nxo)
          call prodel_2(4,wls,mpe,iws,iwa,nxo)
        end do
        !call prodel(4,wls,mpe,iws,iwa)
      end if

      imi = Mul(lsm_inn(lri),ns_sm)
      lrjsta = lri+1
      do lrj=lrjsta,norb_dz
        imj = lsm_inn(lrj)
        imij = Mul(imi,imj)
        if (imij /= imad) cycle
        iws = just(lri,lrj)
        ! s1: d&r&l(1)
        vl0 = fqi*v_onevsqtwo*vlop0
        vl1 = ps1*vlop1
        !wls = -Two*vl0*vdint(lra,lrj)+(vl0-vl1)*voint(lra,lrj)
        wls = -Two*vl0
        call trans_ijkl_intpos(lra,lra,lrj,lrj,nxo)
        call prodel_2(4,wls,mpe,iws,iwa,nxo)
        wls = vl0-vl1
        call trans_ijkl_intpos(lra,lrj,lra,lrj,nxo)
        call prodel_2(4,wls,mpe,iws,iwa,nxo)

        ! s4: d&r&l(2)+c"(1)
        vl0 = fqi*v_onevsqtwo*vlop0
        vl1 = ps4*vlop1
        !wls = wls-Two*vl0*vdint(lra,lri)+(vl0-vl1)*voint(lra,lri)
        wls = -Two*vl0
        call trans_ijkl_intpos(lra,lra,lri,lri,nxo)
        call prodel_2(4,wls,mpe,iws,iwa,nxo)
        wls = vl0-vl1
        call trans_ijkl_intpos(lra,lri,lra,lri,nxo)
        call prodel_2(4,wls,mpe,iws,iwa,nxo)

        vl0 = fqi*v_sqtwo*vlop0
        ! s: d&r&l(3)+c"(2)+c"(1)
        do lr=1,lri-1
          !wls = wls+vl0*vdint(lr,lra)
          wls = vl0
          call trans_ijkl_intpos(lra,lr,lra,lr,nxo)
          call prodel_2(4,wls,mpe,iws,iwa,nxo)
          wls = -Two*vl0
          call trans_ijkl_intpos(lra,lra,lr,lr,nxo)
          call prodel_2(4,wls,mpe,iws,iwa,nxo)
        end do
        ! s: d&r&l(3)+c"(1)
        do lr=lri+1,lrj-1
          !wls = wls+vl0*vdint(lr,lra)
          wls = vl0
          call trans_ijkl_intpos(lra,lr,lra,lr,nxo)
          call prodel_2(4,wls,mpe,iws,iwa,nxo)
          wls = -Two*vl0
          call trans_ijkl_intpos(lra,lra,lr,lr,nxo)
          call prodel_2(4,wls,mpe,iws,iwa,nxo)
        end do
        ! s: d&r&l(3)
        do lr=lrj+1,norb_dz
          !wls = wls+vl0*vdint(lr,lra)
          wls = vl0
          call trans_ijkl_intpos(lra,lr,lra,lr,nxo)
          call prodel_2(4,wls,mpe,iws,iwa,nxo)
          wls = -Two*vl0
          call trans_ijkl_intpos(lra,lra,lr,lr,nxo)
          call prodel_2(4,wls,mpe,iws,iwa,nxo)
        end do
        !call prodel(4,wls,mpe,iws,iwa)
      end do
    end do

    if (jb_sys == 0) return      !any diffrence when jb_sys=1 and jb_sy
    fqi = fg
    do lri=norb_frz+1,norb_dz
      imi = Mul(lsm_inn(lri),ns_sm)
      do lrj=lri+1,norb_dz
        imj = lsm_inn(lrj)
        imij = Mul(imi,imj)
        if (imij /= imad) cycle
        !iws = iws+1
        iws = just(lrj,lri)
        ! s1: d&r&l(1)-c"(2)
        vl0 = fqi*v_onevsqtwo*vlop0
        vl1 = ps3*vlop1
        !wls = -Two*vl0*vdint(lra,lri)+(vl0-vl1)*voint(lra,lri)
        wls = -Two*vl0
        call trans_ijkl_intpos(lra,lra,lri,lri,nxo)
        call prodel_2(4,wls,mpe,iws,iwa,nxo)
        wls = vl0-vl1
        call trans_ijkl_intpos(lra,lri,lra,lri,nxo)
        call prodel_2(4,wls,mpe,iws,iwa,nxo)

        ! s3: (11)d&r&l(2)
        vl0 = fqi*v_onevsqtwo*vlop0
        vl1 = ps2*vlop1
        !wls = wls-Two*vl0*vdint(lra,lrj)+(vl0-vl1)*voint(lra,lrj)
        wls = -Two*vl0
        call trans_ijkl_intpos(lra,lra,lrj,lrj,nxo)
        call prodel_2(4,wls,mpe,iws,iwa,nxo)
        wls = vl0-vl1
        call trans_ijkl_intpos(lra,lrj,lra,lrj,nxo)
        call prodel_2(4,wls,mpe,iws,iwa,nxo)

        vl0 = fqi*v_sqtwo*vlop0
        ! s: d&r&l(3)+c"(1)+c"(2)
        do lr=1,lri-1
          !wls = wls+vl0*vdint(lr,lra)
          wls = vl0
          call trans_ijkl_intpos(lra,lr,lra,lr,nxo)
          call prodel_2(4,wls,mpe,iws,iwa,nxo)
          wls = -Two*vl0
          call trans_ijkl_intpos(lra,lra,lr,lr,nxo)
          call prodel_2(4,wls,mpe,iws,iwa,nxo)
        end do
        ! s: d&r&l(3)+c"(2)
        do lr=lri+1,lrj-1
          !wls = wls+vl0*vdint(lr,lra)
          wls = vl0
          call trans_ijkl_intpos(lra,lr,lra,lr,nxo)
          call prodel_2(4,wls,mpe,iws,iwa,nxo)
          wls = -Two*vl0
          call trans_ijkl_intpos(lra,lra,lr,lr,nxo)
          call prodel_2(4,wls,mpe,iws,iwa,nxo)
        end do
        ! s: d&r&l(3)
        do lr=lrj+1,norb_dz
          !wls = wls+vl0*vdint(lr,lra)
          wls = vl0
          call trans_ijkl_intpos(lra,lr,lra,lr,nxo)
          call prodel_2(4,wls,mpe,iws,iwa,nxo)
          wls = -Two*vl0
          call trans_ijkl_intpos(lra,lra,lr,lr,nxo)
          call prodel_2(4,wls,mpe,iws,iwa,nxo)
        end do
        !call prodel(4,wls,mpe,iws,iwa)
      end do
    end do

  case (5)
    fqi = -fg
    iwd = 0
    do lri=norb_frz+1,norb_dz
      imd = Mul(lsm_inn(lri),ns_sm)
      if (imd /= imad) cycle
      iwd = jud(lri)

      ! dd1: d&r&l(1)
      vl0 = fqi*v_onevsqtwo*vlop0
      vl1 = pdd*vlop1
      !wld = -Two*vl0*vdint(lra,lri)+(vl0-vl1)*voint(lra,lri)
      wld = vl0-vl1
      call trans_ijkl_intpos(lra,lri,lra,lri,nxo)
      call prodel_2(4,wld,mpe,iwd,iwa,nxo)
      wld = -Two*vl0
      call trans_ijkl_intpos(lra,lra,lri,lri,nxo)
      call prodel_2(4,wld,mpe,iwd,iwa,nxo)

      vl0 = fqi*v_sqtwo*vlop0
      ! d: d&r&l(3)+c"(1)
      do lr=1,lri-1
        !wld = wld+vl0*vdint(lr,lra)
        wld = vl0
        call trans_ijkl_intpos(lra,lr,lra,lr,nxo)
        call prodel_2(4,wld,mpe,iwd,iwa,nxo)
        wld = -Two*vl0
        call trans_ijkl_intpos(lra,lra,lr,lr,nxo)
        call prodel_2(4,wld,mpe,iwd,iwa,nxo)
      end do
      ! d: d&r&l(3)
      do lr=lri+1,norb_dz
        !wld = wld+vl0*vdint(lr,lra)
        wld = vl0
        call trans_ijkl_intpos(lra,lr,lra,lr,nxo)
        call prodel_2(4,wld,mpe,iwd,iwa,nxo)
        wld = -Two*vl0
        call trans_ijkl_intpos(lra,lra,lr,lr,nxo)
        call prodel_2(4,wld,mpe,iwd,iwa,nxo)
      end do

      !call prodel(4,wld,mpe,iwd,iwa)
    end do

  case (6)
    fqi = fg
    iwt = 0
    do lri=norb_frz+1,norb_dz
      imi = Mul(lsm_inn(lri),ns_sm)
      do lrj=lri+1,norb_dz
        imj = lsm_inn(lrj)
        imij = Mul(imi,imj)
        if (imij /= imad) cycle
        iwt = just(lri,lrj)

        ! tt: d&r&l(1)
        ! tt: d&r&l(1)+c"(1)
        vl0 = fqi*v_onevsqtwo*vlop0
        vl1 = ptt*vlop1
        !wlt = -Two*vl0*vdint(lra,lri)+(vl0-vl1)*voint(lra,lri)
        wlt = -Two*vl0
        call trans_ijkl_intpos(lra,lra,lri,lri,nxo)
        call prodel_2(4,wlt,mpe,iwt,iwa,nxo)
        wlt = vl0-vl1
        call trans_ijkl_intpos(lra,lri,lra,lri,nxo)
        call prodel_2(4,wlt,mpe,iwt,iwa,nxo)

        !wlt = wlt-Two*vl0*vdint(lra,lrj)+(vl0-vl1)*voint(lra,lrj)
        wlt = -Two*vl0
        call trans_ijkl_intpos(lra,lra,lrj,lrj,nxo)
        call prodel_2(4,wlt,mpe,iwt,iwa,nxo)
        wlt = vl0-vl1
        call trans_ijkl_intpos(lra,lrj,lra,lrj,nxo)
        call prodel_2(4,wlt,mpe,iwt,iwa,nxo)

        vl0 = fqi*v_sqtwo*vlop0
        do lr=1,lri-1
          ! t: d&r&l(3)+c"(1)+c"(1)
          !wlt = wlt+vl0*vdint(lr,lra)
          wlt = vl0
          call trans_ijkl_intpos(lra,lr,lra,lr,nxo)
          call prodel_2(4,wlt,mpe,iwt,iwa,nxo)
          wlt = -Two*vl0
          call trans_ijkl_intpos(lra,lra,lr,lr,nxo)
          call prodel_2(4,wlt,mpe,iwt,iwa,nxo)
        end do
        ! t: d&r&l(3)+c"(1)
        do lr=lri+1,lrj-1
          !wlt = wlt+vl0*vdint(lr,lra)
          wlt = vl0
          call trans_ijkl_intpos(lra,lr,lra,lr,nxo)
          call prodel_2(4,wlt,mpe,iwt,iwa,nxo)
          wlt = -Two*vl0
          call trans_ijkl_intpos(lra,lra,lr,lr,nxo)
          call prodel_2(4,wlt,mpe,iwt,iwa,nxo)
        end do
        ! t: d&r&l(3)
        do lr=lrj+1,norb_dz
          !wlt = wlt+vl0*vdint(lr,lra)
          wlt = vl0
          call trans_ijkl_intpos(lra,lr,lra,lr,nxo)
          call prodel_2(4,wlt,mpe,iwt,iwa,nxo)
          wlt = -Two*vl0
          call trans_ijkl_intpos(lra,lra,lr,lr,nxo)
          call prodel_2(4,wlt,mpe,iwt,iwa,nxo)
        end do

        !call prodel(4,wlt,mpe,iwt,iwa)
      end do
    end do
end select

return

end subroutine diagonal_link_ad_g

subroutine diagonal_link_dae_g(mh)

use gugaci_global, only: fg, jb_sys, jpad, jud, just, jwh, lsm_inn, norb_dz, norb_frz, ns_sm, pd, pdd, ps1, ps2, ps3, ps4, pt, &
                         ptt, th, thh, v_onevsqtwo, v_sqtwo
use Symmetry_Info, only: Mul
use Constants, only: Zero
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: mh
integer(kind=iwp) :: imad, imd, imi, imij, imj, ip, ityad, iwa, iwd, iws, iwt, lri, lrj, lrjsta
real(kind=wp) :: fqi, vij0, vij1, vij2, vl0, vlop0, vlop1

ityad = 1
if (jpad /= 1) ityad = (jpad-1)/8+2
imad = mod(jpad-1,8)
if (imad == 0) then
  ityad = ityad-1
  imad = 8
end if
do ip=1,mh
  iwa = jwh(ip)
  vlop0 = th(ip)
  vlop1 = thh(ip)

  select case (ityad)
    case default ! (1)
      !jpad = 1
      if (abs(vlop0) < 1.0e-30_wp) cycle
      fqi = fg
      ! v: d&r&l(3)
      lri = 0
      lrj = 0
      iwd = 0
      vij0 = Zero
      vij1 = Zero
      vij2 = Zero
      vl0 = fqi*v_sqtwo*vlop0
      call diagonal_call_dae_g(lri,lrj,iwd,iwa,vij0,vij1,vij2,vl0)

    case (2)
      !jpad = jd(im)
      fqi = -fg
      lrj = 0
      do lri=norb_frz+1,norb_dz
        imd = Mul(lsm_inn(lri),ns_sm)
        if (imd /= imad) cycle
        iwd = jud(lri)

        ! d: d&r&l(2)
        vij0 = fqi*v_onevsqtwo*vlop0
        vij1 = pd*vlop1
        vij2 = Zero
        ! d: d&r&l(3)+c"(2)
        ! d: d&r&l(3)
        vl0 = fqi*v_sqtwo*vlop0

        call diagonal_call_dae_g(lri,lrj,iwd,iwa,vij0,vij1,vij2,vl0)

      end do

    case (3)
      !jpad = jt(im)
      fqi = fg
      iwt = 0
      do lri=norb_frz+1,norb_dz
        imi = Mul(lsm_inn(lri),ns_sm)
        do lrj=lri+1,norb_dz
          imj = lsm_inn(lrj)
          imij = Mul(imi,imj)
          if (imij /= imad) cycle
          iwt = just(lri,lrj)

          ! t: d&r&l(2)
          ! t: d&r&l(2)+c"(2)
          vij0 = fqi*v_onevsqtwo*vlop0
          vij1 = pt*vlop1
          vij2 = vij1
          ! t: d&r&l(3)+c"(2)+c"(2)
          ! t: d&r&l(3)+c"(2)
          ! t: d&r&l(3)
          vl0 = fqi*v_sqtwo*vlop0

          call diagonal_call_dae_g(lri,lrj,iwt,iwa,vij0,vij1,vij2,vl0)
        end do
      end do

    case (4)
      fqi = fg
      iws = 0
      do lri=norb_frz+1,norb_dz      !cc
        if (imad == ns_sm) then
          lrj = lri
          iws = just(lri,lri)
          vij0 = Zero
          vij1 = Zero
          vij2 = Zero
          ! s: d&r&l(3)+c"(0)
          ! s: d&r&l(3)
          vl0 = fqi*v_sqtwo*vlop0
          call diagonal_call_dae_g(lri,lrj,iws,iwa,vij0,vij1,vij2,vl0)
        end if

        imi = Mul(lsm_inn(lri),ns_sm)
        lrjsta = lri+1
        do lrj=lrjsta,norb_dz
          imj = lsm_inn(lrj)
          imij = Mul(imi,imj)
          if (imij /= imad) cycle
          iws = just(lri,lrj)
          ! s2: d&r&l(2)
          ! s4: d&r&l(2)+c"(1)
          vij0 = fqi*v_onevsqtwo*vlop0
          vij1 = ps1*vlop1
          vij2 = ps4*vlop1
          ! s: d&r&l(3)+c"(2)+c"(1)
          ! s: d&r&l(3)+c"(1)
          ! s: d&r&l(3)
          vl0 = fqi*v_sqtwo*vlop0
          call diagonal_call_dae_g(lri,lrj,iws,iwa,vij0,vij1,vij2,vl0)
        end do
      end do
      if (jb_sys == 0) cycle
      fqi = fg
      do lri=norb_frz+1,norb_dz
        imi = Mul(lsm_inn(lri),ns_sm)
        if (imad == ns_sm) lrjsta = lri
        do lrj=lri+1,norb_dz
          imj = lsm_inn(lrj)
          imij = Mul(imi,imj)
          if (imij /= imad) cycle
          iws = just(lrj,lri)
          ! s1: d&r&l(1)
          ! s3: d&r&l(1)+c"(2)
          vij0 = fqi*v_onevsqtwo*vlop0
          vij1 = ps2*vlop1
          vij2 = ps3*vlop1
          ! s: d&r&l(3)+c"(1)+c"(2)
          ! s: d&r&l(3)+c"(2)
          ! s: d&r&l(3)
          vl0 = fqi*v_sqtwo*vlop0

          call diagonal_call_dae_g(lri,lrj,iws,iwa,vij0,vij1,vij2,vl0)

        end do
      end do

    case (5)
      fqi = -fg
      lrj = 0
      do lri=norb_frz+1,norb_dz
        imd = Mul(lsm_inn(lri),ns_sm)
        if (imd /= imad) cycle
        iwd = jud(lri)

        ! dd1: d&r&l(1)
        vij0 = fqi*v_onevsqtwo*vlop0
        vij1 = pdd*vlop1
        vij2 = Zero
        ! d: d&r&l(3)+c"(1)
        ! d: d&r&l(3)
        vl0 = fqi*v_sqtwo*vlop0

        call diagonal_call_dae_g(lri,lrj,iwd,iwa,vij0,vij1,vij2,vl0)

      end do

    case (6)
      fqi = fg           !aa
      iwt = 0
      do lri=norb_frz+1,norb_dz
        imi = Mul(lsm_inn(lri),ns_sm)
        do lrj=lri+1,norb_dz
          imj = lsm_inn(lrj)
          imij = Mul(imi,imj)
          if (imij /= imad) cycle
          iwt = just(lri,lrj)

          vij0 = fqi*v_onevsqtwo*vlop0
          ! tt: d&r&l(1)
          vij1 = ptt*vlop1
          ! tt: d&r&l(1)+c"(1)
          vij2 = ptt*vlop1
          ! t: d&r&l(3)+c"(1)+c"(1)
          ! t: d&r&l(3)+c"(1)
          ! t: d&r&l(3)
          vl0 = fqi*v_sqtwo*vlop0
          call diagonal_call_dae_g(lri,lrj,iwt,iwa,vij0,vij1,vij2,vl0)
        end do
      end do
  end select

end do

return

end subroutine diagonal_link_dae_g

subroutine diagonal_call_dae_g(lri,lrj,iwd,iwa,vij0,vij1,vij2,vl0)

use gugaci_global, only: ibsm_ext, iesm_ext, ipae, jpad, ng_sm, norb_all, norb_dz, norb_ext, v_onevsqtwo, v_sqthreevsqtwo, v_sqtwo
use Symmetry_Info, only: Mul
use Constants, only: Two
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: lri, lrj, iwd, iwa
real(kind=wp), intent(in) :: vij0, vij1, vij2, vl0
integer(kind=iwp) :: ima, imae, imb, ityae, iwe, la, lb, lbend, lbsta, lr, lra, lrb, nxo
real(kind=wp) :: vlop0, vlop1, wg13, wg14, wl

if (norb_dz == 0) return
if (ipae == 1) return   !could not link
ityae = (ipae-1)/8+1
imae = mod(ipae-1,8)
if (imae == 0) then
  ityae = ityae-1
  imae = 8
end if
iwe = 0
! 520 = <a,j,k,a>:13,14(ss=3),38(tt=2),50(dd=1)
!link arc_d
select case (ityae)
  case default ! (1)
    !zz = '  g50  '
    do la=ibsm_ext(imae),iesm_ext(imae)
      lra = norb_all-la+1
      iwe = iwe+1
      vlop0 = vl0*v_onevsqtwo
      !wl = Zero
      do lr=1,norb_dz
        if (lr == lri) cycle
        if (lr == lrj) cycle
        !wl = wl+vlop0*vdint(lr,lra)    !db space drl(33)- ext space
        wl = vlop0
        call trans_ijkl_intpos(lra,lr,lra,lr,nxo)
        call prodel_2(6,wl,iwd,iwa,iwe,nxo)
        wl = -Two*vlop0
        call trans_ijkl_intpos(lra,lra,lr,lr,nxo)
        call prodel_2(6,wl,iwd,iwa,iwe,nxo)
      end do
      if (lri >= lrj) then
        vlop0 = vij0*v_onevsqtwo
        vlop1 = -vij1*v_sqthreevsqtwo
        !wl = wl+(vlop0-vlop1)*voint(lra,lri)-Two*vlop0*vdint(lra,lri)
        if (lri /= 0) then
          wl = vlop0-vlop1
          call trans_ijkl_intpos(lra,lri,lra,lri,nxo)
          call prodel_2(6,wl,iwd,iwa,iwe,nxo)
          wl = -Two*vlop0
          call trans_ijkl_intpos(lra,lra,lri,lri,nxo)
          call prodel_2(6,wl,iwd,iwa,iwe,nxo)
        end if
        vlop1 = -vij2*v_sqthreevsqtwo
        !wl = wl+(vlop0-vlop1)*voint(lra,lrj)-Two*vlop0*vdint(lra,lrj)
        if (lrj /= 0) then
          wl = vlop0-vlop1
          call trans_ijkl_intpos(lra,lrj,lra,lrj,nxo)
          call prodel_2(6,wl,iwd,iwa,iwe,nxo)
          wl = -Two*vlop0
          call trans_ijkl_intpos(lra,lra,lrj,lrj,nxo)
          call prodel_2(6,wl,iwd,iwa,iwe,nxo)
        end if
        !write(u6,'(a11,2i3,i6,1x,5f10.4)') zz,lr0,la,jwl,vo(lr0,la),vmd(lr0,la),wg50,wwg50,wl
      else
        vlop0 = vij0*v_onevsqtwo
        vlop1 = -vij1*v_sqthreevsqtwo  !db space (22)drl(11)- ext space -g
        !wl = wl+(vlop0-vlop1)*voint(lra,lrj)-Two*vlop0*vdint(lra,lrj)
        if (lrj /= 0) then
          wl = vlop0-vlop1
          call trans_ijkl_intpos(lra,lrj,lra,lrj,nxo)
          call prodel_2(6,wl,iwd,iwa,iwe,nxo)
          wl = -Two*vlop0
          call trans_ijkl_intpos(lra,lra,lrj,lrj,nxo)
          call prodel_2(6,wl,iwd,iwa,iwe,nxo)
        end if
        vlop1 = -vij2*v_sqthreevsqtwo  !db space drl(22)c"(11)- ext space
        !wl = wl+(vlop0-vlop1)*voint(lra,lri)-Two*vlop0*vdint(lra,lri)
        if (lri /= 0) then
          wl = vlop0-vlop1
          call trans_ijkl_intpos(lra,lri,lra,lri,nxo)
          call prodel_2(6,wl,iwd,iwa,iwe,nxo)
          wl = -Two*vlop0
          call trans_ijkl_intpos(lra,lra,lri,lri,nxo)
          call prodel_2(6,wl,iwd,iwa,iwe,nxo)
        end if
      end if
      !call prodel(6,wl,iwd,iwa,iwe)
    end do

  case (2)
    !zz = '  g38,39  '
    do ima=1,ng_sm
      imb = Mul(ima,imae)
      if (imb > ima) cycle
      do la=ibsm_ext(ima),iesm_ext(ima)
        lra = norb_all-la+1
        lbsta = ibsm_ext(imb)
        lbend = iesm_ext(imb)
        if (ima == imb) lbend = la-1
        do lb=lbsta,lbend
          lrb = norb_all-lb+1
          iwe = iwe+1
          !volalb = Zero
          !vd2lalb = Zero
          vlop0 = -vl0*v_onevsqtwo

          do lr=1,norb_dz
            if (lr == lri) cycle
            if (lr == lrj) cycle
            !volalb = volalb+(voint(lra,lr)+voint(lrb,lr))
            !vd2lalb = vd2lalb-Two*(vdint(lra,lr)+vdint(lrb,lr))
            wl = vlop0
            call trans_ijkl_intpos(lra,lr,lra,lr,nxo)
            call prodel_2(6,wl,iwd,iwa,iwe,nxo)
            call trans_ijkl_intpos(lrb,lr,lrb,lr,nxo)
            call prodel_2(6,wl,iwd,iwa,iwe,nxo)
            wl = -Two*vlop0
            call trans_ijkl_intpos(lra,lra,lr,lr,nxo)
            call prodel_2(6,wl,iwd,iwa,iwe,nxo)
            call trans_ijkl_intpos(lrb,lrb,lr,lr,nxo)
            call prodel_2(6,wl,iwd,iwa,iwe,nxo)

          end do

          !wl = vlop0*(volalb+vd2lalb)
          if (lri >= lrj) then
            vlop0 = -vij0*v_onevsqtwo
            vlop1 = vij1
            !wl = wl+(vlop0-vlop1)*(voint(lra,lri)+voint(lrb,lri))-Two*vlop0*(vdint(lra,lri)+vdint(lrb,lri))
            if (lri /= 0) then
              wl = vlop0-vlop1
              call trans_ijkl_intpos(lra,lri,lra,lri,nxo)
              call prodel_2(6,wl,iwd,iwa,iwe,nxo)
              call trans_ijkl_intpos(lrb,lri,lrb,lri,nxo)
              call prodel_2(6,wl,iwd,iwa,iwe,nxo)
              wl = -Two*vlop0
              call trans_ijkl_intpos(lra,lra,lri,lri,nxo)
              call prodel_2(6,wl,iwd,iwa,iwe,nxo)
              call trans_ijkl_intpos(lrb,lrb,lri,lri,nxo)
              call prodel_2(6,wl,iwd,iwa,iwe,nxo)
            end if
            vlop1 = vij2
            !wl = wl+(vlop0-vlop1)*(voint(lra,lrj)+voint(lrb,lrj))-Two*vlop0*(vdint(lra,lrj)+vdint(lrb,lrj))
            if (lrj /= 0) then
              wl = vlop0-vlop1
              call trans_ijkl_intpos(lra,lrj,lra,lrj,nxo)
              call prodel_2(6,wl,iwd,iwa,iwe,nxo)
              call trans_ijkl_intpos(lrb,lrj,lrb,lrj,nxo)
              call prodel_2(6,wl,iwd,iwa,iwe,nxo)
              wl = -Two*vlop0
              call trans_ijkl_intpos(lra,lra,lrj,lrj,nxo)
              call prodel_2(6,wl,iwd,iwa,iwe,nxo)
              call trans_ijkl_intpos(lrb,lrb,lrj,lrj,nxo)
              call prodel_2(6,wl,iwd,iwa,iwe,nxo)
            end if
          else
            vlop0 = -vij0*v_onevsqtwo
            vlop1 = vij1
            !wl = wl+(vlop0-vlop1)*(voint(lra,lrj)+voint(lrb,lrj))-Two*vlop0*(vdint(lra,lrj)+vdint(lrb,lrj))
            if (lrj /= 0) then
              wl = vlop0-vlop1
              call trans_ijkl_intpos(lra,lrj,lra,lrj,nxo)
              call prodel_2(6,wl,iwd,iwa,iwe,nxo)
              call trans_ijkl_intpos(lrb,lrj,lrb,lrj,nxo)
              call prodel_2(6,wl,iwd,iwa,iwe,nxo)
              wl = -Two*vlop0
              call trans_ijkl_intpos(lra,lra,lrj,lrj,nxo)
              call prodel_2(6,wl,iwd,iwa,iwe,nxo)
              call trans_ijkl_intpos(lrb,lrb,lrj,lrj,nxo)
              call prodel_2(6,wl,iwd,iwa,iwe,nxo)
            end if
            vlop1 = vij2
            !wl = wl+(vlop0-vlop1)*(voint(lra,lri)+voint(lrb,lri))-Two*vlop0*(vdint(lra,lri)+vdint(lrb,lri))
            if (lri /= 0) then
              wl = vlop0-vlop1
              call trans_ijkl_intpos(lra,lri,lra,lri,nxo)
              call prodel_2(6,wl,iwd,iwa,iwe,nxo)
              call trans_ijkl_intpos(lrb,lri,lrb,lri,nxo)
              call prodel_2(6,wl,iwd,iwa,iwe,nxo)
              wl = -Two*vlop0
              call trans_ijkl_intpos(lra,lra,lri,lri,nxo)
              call prodel_2(6,wl,iwd,iwa,iwe,nxo)
              call trans_ijkl_intpos(lrb,lrb,lri,lri,nxo)
              call prodel_2(6,wl,iwd,iwa,iwe,nxo)
            end if
          end if
          !call prodel(6,wl,iwd,iwa,iwe)
        end do
      end do
    end do

  case (3)
    !===========================lyb=====================================
    !zz = '  g14,15  '
    !do ima=1,8
    ! the 8 should be changed to ng_sm
    do ima=1,ng_sm
      imb = Mul(ima,imae)
      if (imb > ima) cycle
      do la=ibsm_ext(ima),iesm_ext(ima)
        lra = norb_all-la+1
        lbsta = ibsm_ext(imb)
        lbend = iesm_ext(imb)
        if (ima == imb) lbend = la-1
        do lb=lbsta,lbend
          lrb = norb_all-lb+1
          iwe = iwe+1
          !volalb = Zero
          !vd2lalb = Zero

          wg14 = -vl0*v_onevsqtwo

          do lr=1,norb_dz
            if (lr == lri) cycle
            if (lr == lrj) cycle
            !volalb = volalb+(voint(lra,lr)+voint(lrb,lr))
            !vd2lalb = vd2lalb-Two*(vdint(lra,lr)+vdint(lrb,lr))
            wl = wg14
            call trans_ijkl_intpos(lra,lr,lra,lr,nxo)
            call prodel_2(6,wl,iwd,iwa,iwe,nxo)
            call trans_ijkl_intpos(lrb,lr,lrb,lr,nxo)
            call prodel_2(6,wl,iwd,iwa,iwe,nxo)
            wl = -Two*wg14
            call trans_ijkl_intpos(lra,lra,lr,lr,nxo)
            call prodel_2(6,wl,iwd,iwa,iwe,nxo)
            call trans_ijkl_intpos(lrb,lrb,lr,lr,nxo)
            call prodel_2(6,wl,iwd,iwa,iwe,nxo)
          end do

          !wl = wg14*(volalb+vd2lalb)

          if ((jpad == 18) .and. (lri == lrj)) cycle
          wg14 = -vij0*v_onevsqtwo
          !wl = wl+wg14*(voint(lra,lri)+voint(lrb,lri))-Two*wg14*(vdint(lra,lri)+vdint(lrb,lri))
          if (lri /= 0) then
            wl = wg14
            call trans_ijkl_intpos(lra,lri,lra,lri,nxo)
            call prodel_2(6,wl,iwd,iwa,iwe,nxo)
            call trans_ijkl_intpos(lrb,lri,lrb,lri,nxo)
            call prodel_2(6,wl,iwd,iwa,iwe,nxo)
            wl = -Two*wg14
            call trans_ijkl_intpos(lra,lra,lri,lri,nxo)
            call prodel_2(6,wl,iwd,iwa,iwe,nxo)
            call trans_ijkl_intpos(lrb,lrb,lri,lri,nxo)
            call prodel_2(6,wl,iwd,iwa,iwe,nxo)
          end if
          !wl = wl+wg14*(voint(lra,lrj)+voint(lrb,lrj))-Two*wg14*(vdint(lra,lrj)+vdint(lrb,lrj))
          if (lrj /= 0) then
            wl = wg14
            call trans_ijkl_intpos(lra,lrj,lra,lrj,nxo)
            call prodel_2(6,wl,iwd,iwa,iwe,nxo)
            call trans_ijkl_intpos(lrb,lrj,lrb,lrj,nxo)
            call prodel_2(6,wl,iwd,iwa,iwe,nxo)
            wl = -Two*wg14
            call trans_ijkl_intpos(lra,lra,lrj,lrj,nxo)
            call prodel_2(6,wl,iwd,iwa,iwe,nxo)
            call trans_ijkl_intpos(lrb,lrb,lrj,lrj,nxo)
            call prodel_2(6,wl,iwd,iwa,iwe,nxo)
          end if
          !call prodel(6,wl,iwd,iwa,iwe)
        end do
      end do
    end do

    if (ipae /= 18) return
    !zz = '  g13     '
    do la=1,norb_ext
      lra = norb_all-la+1

      iwe = iwe+1
      !vovdla = Zero

      wg13 = -vl0*v_sqtwo

      do lr=1,norb_dz
        if (lr == lri) cycle
        if (lr == lrj) cycle
        !vovdla = vovdla+vdint(lr,lra)
        wl = wg13
        call trans_ijkl_intpos(lra,lr,lra,lr,nxo)
        call prodel_2(6,wl,iwd,iwa,iwe,nxo)
        wl = -Two*wg13
        call trans_ijkl_intpos(lra,lra,lr,lr,nxo)
        call prodel_2(6,wl,iwd,iwa,iwe,nxo)
      end do

      !wl = wg13*vovdla

      wg13 = -vij0*v_sqtwo
      !wl = wl+wg13*(vdint(lri,lra)+vdint(lrj,lra))
      if (lri /= 0) then
        wl = wg13
        call trans_ijkl_intpos(lra,lri,lra,lri,nxo)
        call prodel_2(6,wl,iwd,iwa,iwe,nxo)
        wl = -Two*wg13
        call trans_ijkl_intpos(lra,lra,lri,lri,nxo)
        call prodel_2(6,wl,iwd,iwa,iwe,nxo)
      end if
      if (lrj /= 0) then
        wl = wg13
        call trans_ijkl_intpos(lra,lrj,lra,lrj,nxo)
        call prodel_2(6,wl,iwd,iwa,iwe,nxo)
        wl = -Two*wg13
        call trans_ijkl_intpos(lra,lra,lrj,lrj,nxo)
        call prodel_2(6,wl,iwd,iwa,iwe,nxo)
      end if
      !call prodel(6,wl,iwd,iwa,iwe)
    end do
end select

return

end subroutine diagonal_call_dae_g

subroutine diagonal_dbl_g()

use gugaci_global, only: ipae, iw_downwei, jb_sys, jpad, jud, just, lsm_inn, norb_dz, norb_frz, ns_sm, nu_ae
use Symmetry_Info, only: Mul
use Constants, only: Zero, One, Two, Four, Half
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp) :: ipae_, iwa, iwad, iwd, iwdownv, iws, iws1, iwt, jpad1, jpas, jpat, jpat1, lr, lr0, lrd, lrds, lrg, lrm, mr, &
                     mr0, mrm, nxo, nxo1_0, nxo2_0, nxo_1, nxo_2, nxod_1, nxod_2, nxos_1, nxos_2
real(kind=wp) :: db, w1, wld_1, wld_2, wls1_1, wls1_2, wls_1, wls_2, wlt_1, wlt_2, wt0_1, wt0_2
logical(kind=iwp) :: logic_lij
integer(kind=iwp), external :: iwalk_ad

if (norb_dz == 0) return
!wt0 = Zero
do lr=1,norb_dz             ! ........800...
  !wt0 = wt0+voint(lr,lr)+voint(lr,lr)+vdint(lr,lr)
  wt0_1 = Two
  wt0_2 = One
  call trans_ijkl_intpos(lr,lr,lr,lr,nxo)
  do ipae_=1,25
    ipae = ipae_ ! ipae is in global module, is this necessary?
    if (nu_ae(ipae) == 0) cycle
    iwdownv = iw_downwei(1,ipae)
    do iwa=0,iwdownv-1
      !zz = ' doub_800_v'
      iwad = iwalk_ad(1,ipae,iwa,0)
      call prodel_1(1,wt0_1,0,ipae,iwad,lr,lr)
      call prodel_2(1,wt0_2,0,ipae,iwad,nxo)
      !zz = ' doub_800_s'
    end do
  end do
end do

!jps = js(1)

mr0 = 1
do lr0=norb_frz+1,norb_dz
  do lrd=1,norb_dz
    mr0 = Mul(lsm_inn(lr0),ns_sm)
    iwd = jud(lr0)
    jpad = 1+mr0
    jpad1 = jpad+24

    call trans_ijkl_intpos(lrd,lrd,lrd,lrd,nxo)

    if (lr0 == lrd) then
      ! ......d_800.
      !wld = wt0-voint(lr0,lr0)-vdint(lr0,lr0)
      !wls = wld-voint(lr0,lr0)
      wld_1 = One
      do ipae_=1,25
        ipae = ipae_ ! ipae is in global module, is this necessary?
        if (nu_ae(ipae) == 0) cycle
        iwdownv = iw_downwei(jpad,ipae)
        do iwa=0,iwdownv-1
          iwad = iwalk_ad(jpad,ipae,iwa,iwd)
          call prodel_1(1,wld_1,0,ipae,iwad,lrd,lrd)
          !call prodel(1,wld,0,ipae,iwad)
        end do
      end do

      if (jb_sys > 0) then
        do ipae_=1,25
          ipae = ipae_ ! ipae is in global module, is this necessary?
          if (nu_ae(ipae) == 0) cycle
          iwdownv = iw_downwei(jpad1,ipae)
          do iwa=0,iwdownv-1
            iwad = iwalk_ad(jpad1,ipae,iwa,iwd)
            call prodel_1(1,wld_1,0,ipae,iwad,lrd,lrd)
            !call prodel(1,wld,0,ipae,iwad)
          end do
        end do
      end if

    else
      wld_1 = Two
      wld_2 = One
      wls_1 = Two
      wls_2 = One

      do ipae_=1,25
        ipae = ipae_ ! ipae is in global module, is this necessary?
        if (nu_ae(ipae) == 0) cycle
        iwdownv = iw_downwei(jpad,ipae)
        do iwa=0,iwdownv-1
          iwad = iwalk_ad(jpad,ipae,iwa,iwd)
          call prodel_1(1,wld_1,0,ipae,iwad,lrd,lrd)
          call prodel_2(1,wld_2,0,ipae,iwad,nxo)
          !call prodel(1,wld,0,ipae,iwad)
        end do
      end do

      if (jb_sys > 0) then
        do ipae_=1,25
          ipae = ipae_ ! ipae is in global module, is this necessary?
          if (nu_ae(ipae) == 0) cycle
          iwdownv = iw_downwei(jpad1,ipae)
          do iwa=0,iwdownv-1
            iwad = iwalk_ad(jpad1,ipae,iwa,iwd)
            call prodel_1(1,wld_1,0,ipae,iwad,lrd,lrd)
            call prodel_2(1,wld_2,0,ipae,iwad,nxo)
            !call prodel(1,wld,0,ipae,iwad)
          end do
        end do
      end if
      ! 0.....s_800.
      jpad = 17+ns_sm
      iwd = just(lr0,lr0)
      do ipae_=1,25
        ipae = ipae_ ! ipae is in global module, is this necessary?
        if (nu_ae(ipae) == 0) cycle
        iwdownv = iw_downwei(jpad,ipae)
        do iwa=0,iwdownv-1
          iwad = iwalk_ad(jpad,ipae,iwa,iwd)
          call prodel_1(1,wls_1,0,ipae,iwad,lrd,lrd)
          call prodel_2(1,wls_2,0,ipae,iwad,nxo)
          !call prodel(1,wls,0,ipae,iwad)
        end do
      end do
    end if
  end do
  if (lr0 == norb_dz) cycle
  !wld0 = wld
  ! ........800...
  do lr=lr0+1,norb_dz
    mr = Mul(mr0,lsm_inn(lr))
    jpat = 9+mr
    jpas = 17+mr
    jpat1 = jpat+24
    iws = just(lr0,lr)
    iwt = iws
    !wld = wld0-voint(lr,lr)-vdint(lr,lr)
    !=============
    ! lr0 and lr

    wld_1 = One
    wld_2 = Zero
    call trans_ijkl_intpos(lr0,lr0,lr0,lr0,nxo)
    nxo_1 = nxo
    call trans_ijkl_intpos(lr,lr,lr,lr,nxo)
    nxo_2 = nxo

    do ipae_=1,25
      ipae = ipae_ ! ipae is in global module, is this necessary?
      if (nu_ae(ipae) == 0) cycle
      iwdownv = iw_downwei(jpat,ipae)
      do iwa=0,iwdownv-1
        iwad = iwalk_ad(jpat,ipae,iwa,iwt)
        call prodel_1(1,wld_1,0,ipae,iwad,lr0,lr0)
        call prodel_2(1,wld_2,0,ipae,iwad,nxo_1)
        call prodel_1(1,wld_1,0,ipae,iwad,lr,lr)
        call prodel_2(1,wld_2,0,ipae,iwad,nxo_2)
        !call prodel(1,wld,0,ipae,iwad)
      end do
      if (jb_sys > 1) then
        iwdownv = iw_downwei(jpat1,ipae)
        do iwa=0,iwdownv-1
          iwad = iwalk_ad(jpat1,ipae,iwa,iwt)         !t1
          call prodel_1(1,wld_1,0,ipae,iwad,lr0,lr0)
          call prodel_2(1,wld_2,0,ipae,iwad,nxo_1)
          call prodel_1(1,wld_1,0,ipae,iwad,lr,lr)
          call prodel_2(1,wld_2,0,ipae,iwad,nxo_2)
          !call prodel(1,wld,0,ipae,iwad)
        end do
      end if
      iwdownv = iw_downwei(jpas,ipae)
      do iwa=0,iwdownv-1
        iwad = iwalk_ad(jpas,ipae,iwa,iws)
        call prodel_1(1,wld_1,0,ipae,iwad,lr0,lr0)
        call prodel_2(1,wld_2,0,ipae,iwad,nxo_1)
        call prodel_1(1,wld_1,0,ipae,iwad,lr,lr)
        call prodel_2(1,wld_2,0,ipae,iwad,nxo_2)
        !call prodel(1,wld,0,ipae,iwad)
      end do
      if (jb_sys > 0) then
        iws1 = just(lr,lr0)
        do iwa=0,iwdownv-1
          iwad = iwalk_ad(jpas,ipae,iwa,iws1)
          call prodel_1(1,wld_1,0,ipae,iwad,lr0,lr0)
          call prodel_2(1,wld_2,0,ipae,iwad,nxo_1)
          call prodel_1(1,wld_1,0,ipae,iwad,lr,lr)
          call prodel_2(1,wld_2,0,ipae,iwad,nxo_2)
          !call prodel(1,wld,0,ipae,iwad)
        end do
      end if
    end do

    !==============
    ! 1 <= lrd < norb_dz and lrd /= lr0 and lrd /= lr
    do lrd=1,norb_dz

      if ((lrd /= lr0) .and. (lrd /= lr)) then
        wld_1 = Two
        wld_2 = One
        call trans_ijkl_intpos(lrd,lrd,lrd,lrd,nxo)
        do ipae_=1,25
          ipae = ipae_ ! ipae is in global module, is this necessary?
          if (nu_ae(ipae) == 0) cycle
          iwdownv = iw_downwei(jpat,ipae)
          do iwa=0,iwdownv-1
            iwad = iwalk_ad(jpat,ipae,iwa,iwt)
            call prodel_1(1,wld_1,0,ipae,iwad,lrd,lrd)
            call prodel_2(1,wld_2,0,ipae,iwad,nxo)
            !call prodel(1,wld,0,ipae,iwad)
          end do
          if (jb_sys > 1) then
            iwdownv = iw_downwei(jpat1,ipae)
            do iwa=0,iwdownv-1
              iwad = iwalk_ad(jpat1,ipae,iwa,iwt)         !t1
              call prodel_1(1,wld_1,0,ipae,iwad,lrd,lrd)
              call prodel_2(1,wld_2,0,ipae,iwad,nxo)
              !call prodel(1,wld,0,ipae,iwad)
            end do
          end if
          iwdownv = iw_downwei(jpas,ipae)
          do iwa=0,iwdownv-1
            iwad = iwalk_ad(jpas,ipae,iwa,iws)
            call prodel_1(1,wld_1,0,ipae,iwad,lrd,lrd)
            call prodel_2(1,wld_2,0,ipae,iwad,nxo)
            !call prodel(1,wld,0,ipae,iwad)
          end do
          if (jb_sys > 0) then
            iws1 = just(lr,lr0)
            do iwa=0,iwdownv-1
              iwad = iwalk_ad(jpas,ipae,iwa,iws1)
              call prodel_1(1,wld_1,0,ipae,iwad,lrd,lrd)
              call prodel_2(1,wld_2,0,ipae,iwad,nxo)
              !call prodel(1,wld,0,ipae,iwad)
            end do
          end if
        end do
      end if
    end do
  end do
end do

!wl8 = 1/2*hnil*(hnil-1)*vmd(lr,lr)+hnil*voint(lr,lr)  ! 800
!wl5 = (vlop0-vlop1)*vo(lr0,lr)-2*vlop0*vmd(lr0,lr)
!wl5 = vlop0*vmd(lr,lr0)-vlop1*vo(lr0,lr)    ! 2000.11.26
!wt0 = Zero
do lr0=2,norb_dz
  do lr=1,lr0-1
    ! .........520...   vlop0=-2 vlop1=0
    !wt0 = wt0-Two*vdint(lr,lr0)
    wt0_1 = -Two
    wt0_2 = Four
    call trans_ijkl_intpos(lr0,lr,lr0,lr,nxo)
    nxo1_0 = nxo
    call trans_ijkl_intpos(lr0,lr0,lr,lr,nxo)
    nxo2_0 = nxo
    jpad = 1
    iwd = 0
    do ipae_=1,25
      ipae = ipae_ ! ipae is in global module, is this necessary?
      if (nu_ae(ipae) == 0) cycle
      iwdownv = iw_downwei(jpad,ipae)
      do iwa=0,iwdownv-1
        iwad = iwalk_ad(jpad,ipae,iwa,iwd)
        call prodel_2(1,wt0_1,0,ipae,iwad,nxo1_0)
        call prodel_2(1,wt0_2,0,ipae,iwad,nxo2_0)
        !call prodel(1,wt0,0,ipae,iwad)
      end do
    end do
  end do
end do
do lrm=norb_frz+1,norb_dz
  mrm = Mul(lsm_inn(lrm),ns_sm)
  iws = just(lrm,lrm)
  iwd = jud(lrm)
  jpad = 1+mrm
  jpad1 = jpad+24
  jpas = 17+ns_sm
  do lrd=2,norb_dz
    do lrds=1,lrd-1

      logic_lij = .false.
      call trans_ijkl_intpos(lrd,lrds,lrd,lrds,nxo)
      nxod_1 = nxo
      call trans_ijkl_intpos(lrd,lrd,lrds,lrds,nxo)
      nxod_2 = nxo

      do lr=1,lrm-1
        if ((lr == lrds) .and. (lrm == lrd)) then
          logic_lij = .true.
          !wld = wld+vdint(lr,lrm)
          wld_1 = -One
          wld_2 = Two
          do ipae_=1,25
            ipae = ipae_ ! ipae is in global module, is this necessary?
            if (nu_ae(ipae) == 0) cycle
            iwdownv = iw_downwei(jpad,ipae)
            do iwa=0,iwdownv-1
              iwad = iwalk_ad(jpad,ipae,iwa,iwd)
              call prodel_2(1,wld_1,0,ipae,iwad,nxod_1)
              call prodel_2(1,wld_2,0,ipae,iwad,nxod_2)

              !call prodel(1,wld,0,ipae,iwad)
            end do
            if (jb_sys > 0) then
              iwdownv = iw_downwei(jpad1,ipae)
              do iwa=0,iwdownv-1
                iwad = iwalk_ad(jpad1,ipae,iwa,iwd)
                call prodel_2(1,wld_1,0,ipae,iwad,nxod_1)
                call prodel_2(1,wld_2,0,ipae,iwad,nxod_2)

                !call prodel(1,wld,0,ipae,iwad)
              end do
            end if
          end do
          exit
        end if
      end do

      do lr0=lrm+1,norb_dz
        if ((lrm == lrds) .and. (lr0 == lrd)) then
          logic_lij = .true.
          !wld = wld+vdint(lrm,lr0)
          wld_1 = -One
          wld_2 = Two
          do ipae_=1,25
            ipae = ipae_ ! ipae is in global module, is this necessary?
            if (nu_ae(ipae) == 0) cycle
            iwdownv = iw_downwei(jpad,ipae)
            do iwa=0,iwdownv-1
              iwad = iwalk_ad(jpad,ipae,iwa,iwd)
              call prodel_2(1,wld_1,0,ipae,iwad,nxod_1)
              call prodel_2(1,wld_2,0,ipae,iwad,nxod_2)

              !call prodel(1,wld,0,ipae,iwad)
            end do
            if (jb_sys > 0) then
              iwdownv = iw_downwei(jpad1,ipae)
              do iwa=0,iwdownv-1
                iwad = iwalk_ad(jpad1,ipae,iwa,iwd)
                call prodel_2(1,wld_1,0,ipae,iwad,nxod_1)
                call prodel_2(1,wld_2,0,ipae,iwad,nxod_2)

                !call prodel(1,wld,0,ipae,iwad)
              end do
            end if
          end do
          exit
        end if
      end do

      !wls = wld
      !do lr=lrm+1,norb_dz
      !  wls = wls+vdint(lrm,lr)
      !end do
      !do lr0=1,lrm-1
      !  wls = wls+vdint(lr0,lrm)
      !end do

      if (.not. logic_lij) then
        wt0_1 = -Two
        wt0_2 = Four
        do ipae_=1,25
          ipae = ipae_ ! ipae is in global module, is this necessary?
          if (nu_ae(ipae) == 0) cycle
          iwdownv = iw_downwei(jpad,ipae)
          do iwa=0,iwdownv-1
            iwad = iwalk_ad(jpad,ipae,iwa,iwd)
            call prodel_2(1,wt0_1,0,ipae,iwad,nxod_1)
            call prodel_2(1,wt0_2,0,ipae,iwad,nxod_2)

            !call prodel(1,wld,0,ipae,iwad)
          end do
          if (jb_sys > 0) then
            iwdownv = iw_downwei(jpad1,ipae)
            do iwa=0,iwdownv-1
              iwad = iwalk_ad(jpad1,ipae,iwa,iwd)
              call prodel_2(1,wt0_1,0,ipae,iwad,nxod_1)
              call prodel_2(1,wt0_2,0,ipae,iwad,nxod_2)

              !call prodel(1,wld,0,ipae,iwad)
            end do
          end if
          ! ....s_520.
          iwdownv = iw_downwei(jpas,ipae)
          do iwa=0,iwdownv-1
            iwad = iwalk_ad(jpas,ipae,iwa,iws)
            call prodel_2(1,wt0_1,0,ipae,iwad,nxod_1)
            call prodel_2(1,wt0_2,0,ipae,iwad,nxod_2)

            !call prodel(1,wls,0,ipae,iwad)
          end do
        end do
      end if
    end do
  end do

end do

! ........520.
do lr0=norb_frz+1,norb_dz-1
  mr0 = Mul(lsm_inn(lr0),ns_sm)
  do lr=lr0+1,norb_dz
    mr = Mul(mr0,lsm_inn(lr))

    jpat = 9+mr
    jpas = 17+mr
    jpat1 = jpat+24
    iws = just(lr0,lr)
    iwt = iws
    iws1 = just(lr,lr0)

    do lrd=2,norb_dz
      do lrds=1,lrd-1

        logic_lij = .false.
        call trans_ijkl_intpos(lrd,lrds,lrd,lrds,nxo)
        nxos_1 = nxo
        call trans_ijkl_intpos(lrd,lrd,lrds,lrds,nxo)
        nxos_2 = nxo
        if ((lr0 == lrds) .and. (lr == lrd)) then
          logic_lij = .true.
          if (jb_sys == 0) then
            !wls = wt0+Three*(voint(lr,lr0)-vdint(lr,lr0))
            !wlt = wt0+voint(lr,lr0)-Three*vdint(lr,lr0)
            wls_1 = One
            wls_2 = One
            wlt_1 = -One
            wlt_2 = One
          end if
          ! ..3.3  2*vo(lr0,lr)-4*vmd(lr0,lr)  w0=-2  w1=0
          ! ..2.1    vo(lr0,lr)+  vmd(lr0,lr)  w0=-1/2  w1=-3/2
          ! ..3.3  2*vo(lr0,lr)-4*vmd(lr0,lr)  w0=-2  w1=0
          ! ..2.2  - vo(lr0,lr)+  vmd(lr0,lr)  w0=-1/2  w1=1/2
          if (jb_sys > 0) then
            db = jb_sys
            w1 = -(db+3)/(Two*db+2)
            !wls=wt0+(OneHalf-w1)*voint(lr,lr0)-Three*vdint(lr,lr0)
            wls_1 = -Half-w1
            wls_2 = One
            !wlt=wt0+voint(lr,lr0)-Three*vdint(lr,lr0)
            wlt_1 = -One
            wlt_2 = One
            w1 = -(db-1)/(Two*db+2)
            !wls1=wt0+(OneHalf-w1)*voint(lr,lr0)-Three*vdint(lr,lr0)
            wls1_1 = -Half-w1
            wls1_2 = One
            !   ..3.3  2*vo(lr0,lr)-4*vmd(lr0,lr)  w0=-2  w1=0
            !   ..1.2  (-1/2-w1)*vo(lr0,lr)+vmd(lr0,lr)
            ! w0=-1/2    w1=-(db-1)/(2*db+2)
          end if
          do ipae_=1,25
            ipae = ipae_ ! ipae is in global module, is this necessary?
            if (nu_ae(ipae) == 0) cycle
            iwdownv = iw_downwei(jpat,ipae)
            do iwa=0,iwdownv-1
              iwad = iwalk_ad(jpat,ipae,iwa,iwt)
              call prodel_2(1,wlt_1,0,ipae,iwad,nxos_1)
              call prodel_2(1,wlt_2,0,ipae,iwad,nxos_2)
              !call prodel(1,wlt,0,ipae,iwad)
            end do
            iwdownv = iw_downwei(jpas,ipae)
            do iwa=0,iwdownv-1
              iwad = iwalk_ad(jpas,ipae,iwa,iws)
              call prodel_2(1,wls_1,0,ipae,iwad,nxos_1)
              call prodel_2(1,wls_2,0,ipae,iwad,nxos_2)
              !call prodel(1,wls,0,ipae,iwad)
            end do
            if (jb_sys > 0) then
              iwdownv = iw_downwei(jpas,ipae)
              do iwa=0,iwdownv-1
                iwad = iwalk_ad(jpas,ipae,iwa,iws1)     ! s 1-2
                call prodel_2(1,wls1_1,0,ipae,iwad,nxos_1)
                call prodel_2(1,wls1_2,0,ipae,iwad,nxos_2)
                !call prodel(1,wls1,0,ipae,iwad)
              end do
            end if
            if (jb_sys > 1) then
              iwdownv = iw_downwei(jpat1,ipae)
              do iwa=0,iwdownv-1
                iwad = iwalk_ad(jpat1,ipae,iwa,iwt)     !t 1-1
                call prodel_2(1,wlt_1,0,ipae,iwad,nxos_1)
                call prodel_2(1,wlt_2,0,ipae,iwad,nxos_2)
                !call prodel(1,wlt,0,ipae,iwad)
              end do
            end if
          end do
          cycle
        end if

        do lrg=1,lr0-1
          if ((lrg == lrds) .and. (lr0 == lrd)) then
            logic_lij = .true.
            wls_1 = -One
            wls_2 = Two
            wlt_1 = -One
            wlt_2 = Two
            if (jb_sys > 0) then
              wls1_1 = -One
              wls1_2 = Two
            end if
            do ipae_=1,25
              ipae = ipae_ ! ipae is in global module, is this necessary?
              if (nu_ae(ipae) == 0) cycle
              iwdownv = iw_downwei(jpat,ipae)
              do iwa=0,iwdownv-1
                iwad = iwalk_ad(jpat,ipae,iwa,iwt)
                call prodel_2(1,wlt_1,0,ipae,iwad,nxos_1)
                call prodel_2(1,wlt_2,0,ipae,iwad,nxos_2)
                !call prodel(1,wlt,0,ipae,iwad)
              end do
              iwdownv = iw_downwei(jpas,ipae)
              do iwa=0,iwdownv-1
                iwad = iwalk_ad(jpas,ipae,iwa,iws)
                call prodel_2(1,wls_1,0,ipae,iwad,nxos_1)
                call prodel_2(1,wls_2,0,ipae,iwad,nxos_2)
                !call prodel(1,wls,0,ipae,iwad)
              end do
              if (jb_sys > 0) then
                iwdownv = iw_downwei(jpas,ipae)
                do iwa=0,iwdownv-1
                  iwad = iwalk_ad(jpas,ipae,iwa,iws1)     ! s 1-2
                  call prodel_2(1,wls1_1,0,ipae,iwad,nxos_1)
                  call prodel_2(1,wls1_2,0,ipae,iwad,nxos_2)
                  !call prodel(1,wls1,0,ipae,iwad)
                end do
              end if
              if (jb_sys > 1) then
                iwdownv = iw_downwei(jpat1,ipae)
                do iwa=0,iwdownv-1
                  iwad = iwalk_ad(jpat1,ipae,iwa,iwt)     !t 1-1
                  call prodel_2(1,wlt_1,0,ipae,iwad,nxos_1)
                  call prodel_2(1,wlt_2,0,ipae,iwad,nxos_2)
                  !call prodel(1,wlt,0,ipae,iwad)
                end do
              end if
            end do
            exit
          end if

          if ((lrg == lrds) .and. (lr == lrd)) then
            logic_lij = .true.
            wls_1 = -One
            wls_2 = Two
            wlt_1 = -One
            wlt_2 = Two
            if (jb_sys > 0) then
              wls1_1 = -One
              wls1_2 = Two
            end if
            do ipae_=1,25
              ipae = ipae_ ! ipae is in global module, is this necessary?
              if (nu_ae(ipae) == 0) cycle
              iwdownv = iw_downwei(jpat,ipae)
              do iwa=0,iwdownv-1
                iwad = iwalk_ad(jpat,ipae,iwa,iwt)
                call prodel_2(1,wlt_1,0,ipae,iwad,nxos_1)
                call prodel_2(1,wlt_2,0,ipae,iwad,nxos_2)
                !call prodel(1,wlt,0,ipae,iwad)
              end do
              iwdownv = iw_downwei(jpas,ipae)
              do iwa=0,iwdownv-1
                iwad = iwalk_ad(jpas,ipae,iwa,iws)
                call prodel_2(1,wls_1,0,ipae,iwad,nxos_1)
                call prodel_2(1,wls_2,0,ipae,iwad,nxos_2)
                !call prodel(1,wls,0,ipae,iwad)
              end do
              if (jb_sys > 0) then
                iwdownv = iw_downwei(jpas,ipae)
                do iwa=0,iwdownv-1
                  iwad = iwalk_ad(jpas,ipae,iwa,iws1)     ! s 1-2
                  call prodel_2(1,wls1_1,0,ipae,iwad,nxos_1)
                  call prodel_2(1,wls1_2,0,ipae,iwad,nxos_2)
                  !call prodel(1,wls1,0,ipae,iwad)
                end do
              end if
              if (jb_sys > 1) then
                iwdownv = iw_downwei(jpat1,ipae)
                do iwa=0,iwdownv-1
                  iwad = iwalk_ad(jpat1,ipae,iwa,iwt)     !t 1-1
                  call prodel_2(1,wlt_1,0,ipae,iwad,nxos_1)
                  call prodel_2(1,wlt_2,0,ipae,iwad,nxos_2)
                  !call prodel(1,wlt,0,ipae,iwad)
                end do
              end if
            end do
            exit
          end if
        end do

        do lrg=lr0+1,lr-1
          if ((lr0 == lrds) .and. (lrg == lrd)) then
            logic_lij = .true.
            wls_1 = -One
            wls_2 = Two
            wlt_1 = -One
            wlt_2 = Two
            if (jb_sys > 0) then
              wls1_1 = -One
              wls1_2 = Two
            end if
            do ipae_=1,25
              ipae = ipae_ ! ipae is in global module, is this necessary?
              if (nu_ae(ipae) == 0) cycle
              iwdownv = iw_downwei(jpat,ipae)
              do iwa=0,iwdownv-1
                iwad = iwalk_ad(jpat,ipae,iwa,iwt)
                call prodel_2(1,wlt_1,0,ipae,iwad,nxos_1)
                call prodel_2(1,wlt_2,0,ipae,iwad,nxos_2)
                !call prodel(1,wlt,0,ipae,iwad)
              end do
              iwdownv = iw_downwei(jpas,ipae)
              do iwa=0,iwdownv-1
                iwad = iwalk_ad(jpas,ipae,iwa,iws)
                call prodel_2(1,wls_1,0,ipae,iwad,nxos_1)
                call prodel_2(1,wls_2,0,ipae,iwad,nxos_2)
                !call prodel(1,wls,0,ipae,iwad)
              end do
              if (jb_sys > 0) then
                iwdownv = iw_downwei(jpas,ipae)
                do iwa=0,iwdownv-1
                  iwad = iwalk_ad(jpas,ipae,iwa,iws1)     ! s 1-2
                  call prodel_2(1,wls1_1,0,ipae,iwad,nxos_1)
                  call prodel_2(1,wls1_2,0,ipae,iwad,nxos_2)
                  !call prodel(1,wls1,0,ipae,iwad)
                end do
              end if
              if (jb_sys > 1) then
                iwdownv = iw_downwei(jpat1,ipae)
                do iwa=0,iwdownv-1
                  iwad = iwalk_ad(jpat1,ipae,iwa,iwt)     !t 1-1
                  call prodel_2(1,wlt_1,0,ipae,iwad,nxos_1)
                  call prodel_2(1,wlt_2,0,ipae,iwad,nxos_2)
                  !call prodel(1,wlt,0,ipae,iwad)
                end do
              end if
            end do
            exit
          end if

          if ((lrg == lrds) .and. (lr == lrd)) then
            logic_lij = .true.
            wls_1 = -One
            wls_2 = Two
            wlt_1 = -One
            wlt_2 = Two
            if (jb_sys > 0) then
              wls1_1 = -One
              wls1_2 = Two
            end if
            do ipae_=1,25
              ipae = ipae_ ! ipae is in global module, is this necessary?
              if (nu_ae(ipae) == 0) cycle
              iwdownv = iw_downwei(jpat,ipae)
              do iwa=0,iwdownv-1
                iwad = iwalk_ad(jpat,ipae,iwa,iwt)
                call prodel_2(1,wlt_1,0,ipae,iwad,nxos_1)
                call prodel_2(1,wlt_2,0,ipae,iwad,nxos_2)
                !call prodel(1,wlt,0,ipae,iwad)
              end do
              iwdownv = iw_downwei(jpas,ipae)
              do iwa=0,iwdownv-1
                iwad = iwalk_ad(jpas,ipae,iwa,iws)
                call prodel_2(1,wls_1,0,ipae,iwad,nxos_1)
                call prodel_2(1,wls_2,0,ipae,iwad,nxos_2)
                !call prodel(1,wls,0,ipae,iwad)
              end do
              if (jb_sys > 0) then
                iwdownv = iw_downwei(jpas,ipae)
                do iwa=0,iwdownv-1
                  iwad = iwalk_ad(jpas,ipae,iwa,iws1)     ! s 1-2
                  call prodel_2(1,wls1_1,0,ipae,iwad,nxos_1)
                  call prodel_2(1,wls1_2,0,ipae,iwad,nxos_2)
                  !call prodel(1,wls1,0,ipae,iwad)
                end do
              end if
              if (jb_sys > 1) then
                iwdownv = iw_downwei(jpat1,ipae)
                do iwa=0,iwdownv-1
                  iwad = iwalk_ad(jpat1,ipae,iwa,iwt)     !t 1-1
                  call prodel_2(1,wlt_1,0,ipae,iwad,nxos_1)
                  call prodel_2(1,wlt_2,0,ipae,iwad,nxos_2)
                  !call prodel(1,wlt,0,ipae,iwad)
                end do
              end if
            end do
            exit
          end if
        end do

        do lrg=lr+1,norb_dz
          if ((lr == lrds) .and. (lrg == lrd)) then
            logic_lij = .true.
            wls_1 = -One
            wls_2 = Two
            wlt_1 = -One
            wlt_2 = Two
            if (jb_sys > 0) then
              wls1_1 = -One
              wls1_2 = Two
            end if
            do ipae_=1,25
              ipae = ipae_ ! ipae is in global module, is this necessary?
              if (nu_ae(ipae) == 0) cycle
              iwdownv = iw_downwei(jpat,ipae)
              do iwa=0,iwdownv-1
                iwad = iwalk_ad(jpat,ipae,iwa,iwt)
                call prodel_2(1,wlt_1,0,ipae,iwad,nxos_1)
                call prodel_2(1,wlt_2,0,ipae,iwad,nxos_2)
                !call prodel(1,wlt,0,ipae,iwad)
              end do
              iwdownv = iw_downwei(jpas,ipae)
              do iwa=0,iwdownv-1
                iwad = iwalk_ad(jpas,ipae,iwa,iws)
                call prodel_2(1,wls_1,0,ipae,iwad,nxos_1)
                call prodel_2(1,wls_2,0,ipae,iwad,nxos_2)
                !call prodel(1,wls,0,ipae,iwad)
              end do
              if (jb_sys > 0) then
                iwdownv = iw_downwei(jpas,ipae)
                do iwa=0,iwdownv-1
                  iwad = iwalk_ad(jpas,ipae,iwa,iws1)     ! s 1-2
                  call prodel_2(1,wls1_1,0,ipae,iwad,nxos_1)
                  call prodel_2(1,wls1_2,0,ipae,iwad,nxos_2)
                  !call prodel(1,wls1,0,ipae,iwad)
                end do
              end if
              if (jb_sys > 1) then
                iwdownv = iw_downwei(jpat1,ipae)
                do iwa=0,iwdownv-1
                  iwad = iwalk_ad(jpat1,ipae,iwa,iwt)     !t 1-1
                  call prodel_2(1,wlt_1,0,ipae,iwad,nxos_1)
                  call prodel_2(1,wlt_2,0,ipae,iwad,nxos_2)
                  !call prodel(1,wlt,0,ipae,iwad)
                end do
              end if
            end do
            exit
          end if

          if ((lr0 == lrds) .and. (lrg == lrd)) then
            logic_lij = .true.
            wls_1 = -One
            wls_2 = Two
            wlt_1 = -One
            wlt_2 = Two
            if (jb_sys > 0) then
              wls1_1 = -One
              wls1_2 = Two
            end if
            do ipae_=1,25
              ipae = ipae_ ! ipae is in global module, is this necessary?
              if (nu_ae(ipae) == 0) cycle
              iwdownv = iw_downwei(jpat,ipae)
              do iwa=0,iwdownv-1
                iwad = iwalk_ad(jpat,ipae,iwa,iwt)
                call prodel_2(1,wlt_1,0,ipae,iwad,nxos_1)
                call prodel_2(1,wlt_2,0,ipae,iwad,nxos_2)
                !call prodel(1,wlt,0,ipae,iwad)
              end do
              iwdownv = iw_downwei(jpas,ipae)
              do iwa=0,iwdownv-1
                iwad = iwalk_ad(jpas,ipae,iwa,iws)
                call prodel_2(1,wls_1,0,ipae,iwad,nxos_1)
                call prodel_2(1,wls_2,0,ipae,iwad,nxos_2)
                !call prodel(1,wls,0,ipae,iwad)
              end do
              if (jb_sys > 0) then
                iwdownv = iw_downwei(jpas,ipae)
                do iwa=0,iwdownv-1
                  iwad = iwalk_ad(jpas,ipae,iwa,iws1)     ! s 1-2
                  call prodel_2(1,wls1_1,0,ipae,iwad,nxos_1)
                  call prodel_2(1,wls1_2,0,ipae,iwad,nxos_2)
                  !call prodel(1,wls1,0,ipae,iwad)
                end do
              end if
              if (jb_sys > 1) then
                iwdownv = iw_downwei(jpat1,ipae)
                do iwa=0,iwdownv-1
                  iwad = iwalk_ad(jpat1,ipae,iwa,iwt)     !t 1-1
                  call prodel_2(1,wlt_1,0,ipae,iwad,nxos_1)
                  call prodel_2(1,wlt_2,0,ipae,iwad,nxos_2)
                  !call prodel(1,wlt,0,ipae,iwad)
                end do
              end if
            end do
            exit
          end if
        end do

        if (.not. logic_lij) then
          if (jb_sys == 0) then
            wls_1 = -Two
            wls_2 = Four
            wlt_1 = -Two
            wlt_2 = Four
          end if

          if (jb_sys > 0) then
            wls_1 = -Two
            wls_2 = Four
            wlt_1 = -Two
            wlt_2 = Four
            wls1_1 = -Two
            wls1_2 = Four
          end if

          do ipae_=1,25
            ipae = ipae_ ! ipae is in global module, is this necessary?
            if (nu_ae(ipae) == 0) cycle
            iwdownv = iw_downwei(jpat,ipae)
            do iwa=0,iwdownv-1
              iwad = iwalk_ad(jpat,ipae,iwa,iwt)
              call prodel_2(1,wlt_1,0,ipae,iwad,nxos_1)
              call prodel_2(1,wlt_2,0,ipae,iwad,nxos_2)
              !call prodel(1,wlt,0,ipae,iwad)
            end do
            iwdownv = iw_downwei(jpas,ipae)
            do iwa=0,iwdownv-1
              iwad = iwalk_ad(jpas,ipae,iwa,iws)
              call prodel_2(1,wls_1,0,ipae,iwad,nxos_1)
              call prodel_2(1,wls_2,0,ipae,iwad,nxos_2)
              !call prodel(1,wls,0,ipae,iwad)
            end do
            if (jb_sys > 0) then
              iwdownv = iw_downwei(jpas,ipae)
              do iwa=0,iwdownv-1
                iwad = iwalk_ad(jpas,ipae,iwa,iws1)     ! s 1-2
                call prodel_2(1,wls1_1,0,ipae,iwad,nxos_1)
                call prodel_2(1,wls1_2,0,ipae,iwad,nxos_2)
                !call prodel(1,wls1,0,ipae,iwad)
              end do
            end if
            if (jb_sys > 1) then
              iwdownv = iw_downwei(jpat1,ipae)
              do iwa=0,iwdownv-1
                iwad = iwalk_ad(jpat1,ipae,iwa,iwt)     !t 1-1
                call prodel_2(1,wlt_1,0,ipae,iwad,nxos_1)
                call prodel_2(1,wlt_2,0,ipae,iwad,nxos_2)
                !call prodel(1,wlt,0,ipae,iwad)
              end do
            end if
          end do

        end if
      end do
    end do
  end do
end do
! ------------- end of .......h_delm --------------

return

end subroutine diagonal_dbl_g

subroutine diagonal_ext_g()

use gugaci_global, only: ibsm_ext, iesm_ext, ipae, lsm, ng_sm, norb_all, norb_ext
use Symmetry_Info, only: Mul
use Constants, only: One, Two
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp) :: im, ima, imb, ipas, ipat, jw, jweis, jws, jws0, jwt, la, laend, lasta, lb, lra, lrb, lrzz, mr, mra, nxo
real(kind=wp) :: wld, wls, wlt

jws0 = 0
!do mra=1,8
do mra=1,ng_sm
  ipae = 1+mra
  lasta = ibsm_ext(mra)
  laend = iesm_ext(mra)
  lrzz = laend-lasta+1
  lrzz = lrzz*(lrzz-1)/2
  jws0 = jws0+lrzz
  jw = 0
  do la=lasta,laend
    !jpd = jd(mra)
    jw = jw+1
    lra = norb_all-la+1
    !wld = voint(lra,lra)
    wld = One
    call prodel_1(2,wld,0,ipae,jw,lra,lra)
    !call prodel(2,wld,0,ipae,jw)
  end do
end do

!zz = ' out_800_s'
!jps = js(1)

jweis = jws0
do la=1,norb_ext
  !jpd = jd(mra)
  lra = norb_all-la+1
  jweis = jweis+1
  !wls = Two*voint(lra,lra)+vdint(lra,lra)
  !call prodel(2,wls,0,18,jweis)
  wls = Two
  call prodel_1(2,wls,0,18,jweis,lra,lra)
  wls = One
  call trans_ijkl_intpos(lra,lra,lra,lra,nxo)
  call prodel_2(2,wls,0,18,jweis,nxo)
end do
!do im=1,8
do im=1,ng_sm
  jws = 0
  jwt = 0
  ipat = 9+im
  ipas = 17+im
  do la=2,norb_ext
    lra = norb_all-la+1
    ima = lsm(la)
    do lb=1,la-1
      lrb = norb_all-lb+1
      imb = lsm(lb)
      mr = Mul(ima,imb)
      if (mr /= im) cycle
      !jps = js(mr)
      !jpt = jt(mr)
      jws = jws+1
      jwt = jwt+1
      !wls = voint(lra,lra)+voint(lrb,lrb)
      wls = One
      wlt = wls
      call prodel_1(2,wls,0,ipas,jws,lra,lra)
      call prodel_1(2,wlt,0,ipat,jwt,lra,lra)
      call prodel_1(2,wls,0,ipas,jws,lrb,lrb)
      call prodel_1(2,wlt,0,ipat,jwt,lrb,lrb)

      !wls = wls+voint(lrb,lra)+vdint(lrb,lra)   ! w0=-1/2  w1=-3/2
      wls = One
      call trans_ijkl_intpos(lrb,lra,lrb,lra,nxo)
      call prodel_2(2,wls,0,ipas,jws,nxo)
      call trans_ijkl_intpos(lrb,lrb,lra,lra,nxo)
      call prodel_2(2,wls,0,ipas,jws,nxo)
      !wlt = wlt-voint(lrb,lra)+vdint(lrb,lra)   ! w0=-1/2  w1=1/2
      wlt = -One
      call trans_ijkl_intpos(lrb,lra,lrb,lra,nxo)
      call prodel_2(2,wlt,0,ipat,jwt,nxo)
      wlt = One
      call trans_ijkl_intpos(lrb,lrb,lra,lra,nxo)
      call prodel_2(2,wlt,0,ipat,jwt,nxo)

      !call prodel(2,wls,0,ipas,jws)
      !call prodel(2,wlt,0,ipat,jwt)
    end do
  end do
end do

return

end subroutine diagonal_ext_g

! idb=1  in dbl_space         ity_up=0-5              jd_type,jd_im,iwd
! idb=2  in ext_space         ity_down=0-3            je_type,je_im,iwe
! idb=3  in act_space         ity_up=0-5,itdown=0,3   jp,     mpe,  iwa
! idb=4  between dbl and act  ity_up=0-5,itdown=0,3   mpe,    iwd,  iwa
! idb=5  between act and ext  ity_down=0-3            jp,     iwa,  iwe
! idb=6  between dbl and ext  ity_down=0-3            iwd,    iwa,  iwe

! this subroutine prodel_1 does the dm1 part, which corresponds to voin
subroutine prodel_1(idb,wl,mg1,mg2,mg3,mg6,mg7)

use gugaci_global, only: dm1tmp, ican_a, ihy, ipae, iseg_downwei, iw_downwei, iy, jpad, jpad_upwei, jphy, mxnode, nu_ad, vector1
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: idb, mg1, mg2, mg3, mg6, mg7
real(kind=wp), intent(in) :: wl
integer(kind=iwp) :: ii, in_, isegdownwei, iw, iwa, iwa0, iwad, iwd, iwe, iwupwei, jdbl, jp, jph, jw, jwd, jwnu, jwu, lwnu, mg67, &
                     mm, mpe
integer(kind=iwp), external :: iwalk_ad

!ndr = 88
select case (idb)
  case default ! (1)
    ! in dbl_space
    ipae = mg2
    iwad = mg3
    isegdownwei = iseg_downwei(ipae)
    do mm=iwad+1,iwad+isegdownwei
      !vector1(mm) = vector1(mm)+wl
      mg67 = ican_a(mg7)+mg6
      dm1tmp(mg67) = dm1tmp(mg67)+vector1(mm)*wl*vector1(mm)
      !if (mg7 == 4) write(nf2,'(3i4,f18.10)') mg7,mg6,mm,wl
      !write(nf2,'(3i4,3f18.10)') mg7,mg6,mm,wl,vector1(mm),vector1(mm
      !if (mm == ndr) then
      !  write(nf2,'(a8,3i6,2f20.14)') ' in dbl _',mg1,mg2,mg3,wl,vector1(mm)
      !  write(u6,'(a8,3i6,2f20.14)') ' in dbl _',mg1,mg2,mg3,wl,vector1(mm)
      !end if
    end do

  case (2)
    ! in ext_space
    ipae = mg2
    iwe = mg3
    do jdbl=1,mxnode
      if (nu_ad(jdbl) == 0) cycle
      iw = iw_downwei(jdbl,ipae)
      iwupwei = jpad_upwei(jdbl)
      do iwa=0,iw-1
        do iwd=0,iwupwei-1
          mm = iwalk_ad(jdbl,ipae,iwa,iwd)+iwe
          !vector1(mm) = vector1(mm)+wl
          mg67 = ican_a(mg7)+mg6
          dm1tmp(mg67) = dm1tmp(mg67)+vector1(mm)*wl*vector1(mm)

          !if (mm == ndr) then
          !  write(nf2,'(a8,3i6,2f20.14)') ' in ext _',mg1,mg2,mg3,wl,vector1(mm)
          !  write(u6,'(a8,3i6,2f20.14)') ' in ext _',mg1,mg2,mg3,wl,vector1(mm)
          !end if
        end do
      end do
    end do

  case (3)
    ! in act_space
    jp = mg1
    mpe = mg2
    jw = mg3
    iwupwei = jpad_upwei(jpad)
    isegdownwei = iseg_downwei(ipae)
    jph = jphy(jp)
    in_ = ihy(jph)
    lwnu = iy(1,mpe)
    do jwu=jph+1,jph+in_
      iwa = jw+ihy(jwu)-1
      do jwd=1,lwnu
        iwa = iwa+1
        do iwd=0,iwupwei-1
          iwad = iwalk_ad(jpad,ipae,iwa,iwd)
          do iwe=1,isegdownwei
            mm = iwe+iwad
            !vector1(mm) = vector1(mm)+wl
            mg67 = ican_a(mg7)+mg6
            dm1tmp(mg67) = dm1tmp(mg67)+vector1(mm)*wl*vector1(mm)

            !if (mm == ndr) then
            !  write(nf2,'(a8,3i6,2f20.14)') ' in act _',mg1,mg2,mg3,wl,vector1(mm)
            !  write(u6,'(a8,3i6,2f20.14)') ' in act _',mg1,mg2,mg3,wl,vector1(mm)
            !end if
          end do
        end do
      end do
    end do

  case (4)
    ! between dbl and act
    mpe = mg1
    iwd = mg2
    iwa = mg3-1
    isegdownwei = iseg_downwei(ipae)   ! between dbl and act
    jwnu = iy(1,mpe)
    do ii=1,jwnu
      iwa = iwa+1
      iwad = iwalk_ad(jpad,ipae,iwa,iwd)
      do iwe=1,isegdownwei
        mm = iwe+iwad                  ! iwl=iwalk_ad
        !vector1(mm) = vector1(mm)+wl
        mg67 = ican_a(mg7)+mg6
        dm1tmp(mg67) = dm1tmp(mg67)+vector1(mm)*wl*vector1(mm)

        !if (mm == ndr) then
        !  write(nf2,'(a8,3i6,2f20.14)') ' dbl_act ',mg1,mg2,mg3,wl,vector1(mm)
        !  write(u6,'(a8,3i6,2f20.14)') ' dbl_act ',mg1,mg2,mg3,wl,vector1(mm)
        !end if
      end do
    end do

  case (5)
    ! between act and ext
    jp = mg1
    iwa0 = mg2
    iwe = mg3
    iwupwei = jpad_upwei(jpad)
    jph = jphy(jp)
    in_ = ihy(jph)
    do jwu=jph+1,jph+in_
      iwa = iwa0+ihy(jwu)
      do iwd=0,iwupwei-1
        iwad = iwalk_ad(jpad,ipae,iwa,iwd)
        mm = iwe+iwad
        !vector1(mm) = vector1(mm)+wl
        mg67 = ican_a(mg7)+mg6
        dm1tmp(mg67) = dm1tmp(mg67)+vector1(mm)*wl*vector1(mm)

        !if (mm == ndr) then
        !  write(nf2,'(a8,3i6,2f20.14)') ' act_ext ',mg1,mg2,mg3,wl,vector1(mm)
        !  write(u6,'(a8,3i6,2f20.14)') ' act_ext ',mg1,mg2,mg3,wl,vector1(mm)
        !end if
      end do
    end do

  case (6)
    ! between dbl and ext
    iwd = mg1
    iwa = mg2
    iwe = mg3
    iwad = iwalk_ad(jpad,ipae,iwa,iwd)   ! between dbl,act and ext
    mm = iwe+iwad
    !vector1(mm) = vector1(mm)+wl
    mg67 = ican_a(mg7)+mg6
    dm1tmp(mg67) = dm1tmp(mg67)+vector1(mm)*wl*vector1(mm)

    !if (mm == ndr) then
    !  write(nf2,'(a8,3i6,2f20.14)') ' dbl_ext ',mg1,mg2,mg3,wl,vector1(mm)
    !  write(u6,'(a8,3i6,2f20.14)') ' dbl_ext ',mg1,mg2,mg3,wl,vector1(mm)
    !end if
end select

return

end subroutine prodel_1

!this subroutine prodel_2 does the dm2 part, which corresponds to vint_c
subroutine prodel_2(idb,wl,mg1,mg2,mg3,mg6)

use gugaci_global, only: ihy, ipae, iseg_downwei, iw_downwei, iy, jpad, jpad_upwei, jphy, mxnode, nu_ad, vector1, vector2
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: idb, mg1, mg2, mg3, mg6
real(kind=wp), intent(in) :: wl
integer(kind=iwp) :: ii, in_, isegdownwei, iw, iwa, iwa0, iwad, iwd, iwe, iwupwei, jdbl, jp, jph, jw, jwd, jwnu, jwu, lwnu, mm, mpe
integer(kind=iwp), external :: iwalk_ad

!ndr = 88
!ndr = 6
select case (idb)
  case default ! (1)
    ! in dbl_space
    ipae = mg2
    iwad = mg3
    isegdownwei = iseg_downwei(ipae)
    do mm=iwad+1,iwad+isegdownwei
      !vector1(mm) = vector1(mm)+wl
      vector2(mg6) = vector2(mg6)+vector1(mm)*wl*vector1(mm)
      !if (mg6 == 45) write(nf2,'(i8,2i4,4f18.10)') mg6,mm,mm,vector2(mg6),vector1(mm),wl,vector1(mm)

      !if (mg6 == ndr) then
      !  write(nf2,'(a9,3i6,3f20.14)') ' in dbl_ ',mg6,mm,mm,wl,vector1(mm),vector2(mg6)
      !  write(u6,'(a8,3i6,2f20.14)') ' in dbl _',mg1,mg2,mg3,wl,vector1(mm)
      !end if
    end do

  case (2)
    ! in ext_space
    ipae = mg2
    iwe = mg3
    do jdbl=1,mxnode
      if (nu_ad(jdbl) == 0) cycle
      iw = iw_downwei(jdbl,ipae)
      iwupwei = jpad_upwei(jdbl)
      do iwa=0,iw-1
        do iwd=0,iwupwei-1
          mm = iwalk_ad(jdbl,ipae,iwa,iwd)+iwe
          !vector1(mm)=vector1(mm)+wl
          vector2(mg6) = vector2(mg6)+vector1(mm)*wl*vector1(mm)
          !write(nf2,'(i8,2i4,4f18.10)') mg6,mm,mm,vector2(mg6),vector1(mm),wl,vector1(mm)

          !if (mg6 == 2926) then
          !  write(nf2,'(a9,3i6,3f20.14)') ' in ext_ ',mg6,mm,mm,wl,vector1(mm),vector2(mg6)
          !  write(u6,'(a8,3i6,2f20.14)') ' in ext _',mg1,mg2,mg3,wl,vector1(mm)
          !end if
        end do
      end do
    end do

  case (3)
    ! in act_space
    jp = mg1
    mpe = mg2
    jw = mg3
    iwupwei = jpad_upwei(jpad)
    isegdownwei = iseg_downwei(ipae)
    jph = jphy(jp)
    in_ = ihy(jph)
    lwnu = iy(1,mpe)
    do jwu=jph+1,jph+in_
      iwa = jw+ihy(jwu)-1
      do jwd=1,lwnu
        iwa = iwa+1
        do iwd=0,iwupwei-1
          iwad = iwalk_ad(jpad,ipae,iwa,iwd)
          do iwe=1,isegdownwei
            mm = iwe+iwad
            !vector1(mm) = vector1(mm)+wl
            vector2(mg6) = vector2(mg6)+vector1(mm)*wl*vector1(mm)
            !if (mg6 == 231)  write(nf2,'(i4,3f18.10)') mm,wl,vector1(mm),vector2(mg6)
            !write(nf2,'(a3,i4)') 'act',mm

            !write(nf2,'(i8,2i4,4f18.10)') mg6,mm,mm,vector2(mg6),vector1(mm),wl,vector1(mm)

            !if (mg6 == ndr) then
            !  write(nf2,'(a9,3i6,3f20.14)') ' in act_ ',mg6,mm,mm,wl,vector1(mm),vector2(mg6)
            !  write(u6,'(a8,3i6,2f20.14)') ' in act _',mg1,mg2,mg3,wl,vector1(mm)
            !end if
          end do
        end do
      end do
    end do

  case (4)
    ! between dbl and act
    mpe = mg1
    iwd = mg2
    iwa = mg3-1
    isegdownwei = iseg_downwei(ipae)   ! between dbl and act
    jwnu = iy(1,mpe)
    do ii=1,jwnu
      iwa = iwa+1
      iwad = iwalk_ad(jpad,ipae,iwa,iwd)
      do iwe=1,isegdownwei
        mm = iwe+iwad                  ! iwl=iwalk_ad
        !vector1(mm) = vector1(mm)+wl
        vector2(mg6) = vector2(mg6)+vector1(mm)*wl*vector1(mm)
        !if (mg6 == 91) write(nf2,'(i8,i4,3f18.10)') mg6,mm,vector2(mg6),wl,vector1(mm)
        !write(nf2,'(a7,i4)') 'dbl_act',mm

        !if (mg6 == 3) write(nf2,'(i4,2f18.10)') mm,wl,vector1(mm)
        !if (mg6 == ndr) then
        !  write(nf2,'(a9,3i6,3f20.14)') ' dbl_act ',mg6,mm,mm,wl,vector1(mm),vector2(mg6)
        !  write(u6,'(a8,3i6,2f20.14)') ' dbl_act ',mg1,mg2,mg3,wl,vector1(mm)
        !end if
      end do
    end do

  case (5)
    ! between act and ext
    jp = mg1
    iwa0 = mg2
    iwe = mg3
    iwupwei = jpad_upwei(jpad)
    jph = jphy(jp)
    in_ = ihy(jph)
    do jwu=jph+1,jph+in_
      iwa = iwa0+ihy(jwu)
      do iwd=0,iwupwei-1
        iwad = iwalk_ad(jpad,ipae,iwa,iwd)
        mm = iwe+iwad
        !vector1(mm) = vector1(mm)+wl
        vector2(mg6) = vector2(mg6)+vector1(mm)*wl*vector1(mm)
        !write(nf2,*) 'ae',mm,wl,vector1(mm)
        !write(nf2,'(i8,2i4,4f18.10)') mg6,mm,mm,vector2(mg6),vector1(mm),wl,vector1(mm)

        !if (mg6 == ndr) then
        !  write(nf2,'(a9,3i6,3f20.14)') ' act_ext ',mg6,mm,mm,wl,vector1(mm),vector2(mg6)
        !  write(u6,'(a8,3i6,2f20.14)') ' act_ext ',mg1,mg2,mg3,wl,vector1(mm)
        !end if
      end do
    end do

  case (6)
    ! between dbl and ext
    iwd = mg1
    iwa = mg2
    iwe = mg3
    iwad = iwalk_ad(jpad,ipae,iwa,iwd)   ! between dbl,act and ext
    mm = iwe+iwad
    !vector1(mm) = vector1(mm)+wl
    vector2(mg6) = vector2(mg6)+vector1(mm)*wl*vector1(mm)
    !write(nf2,'(i8,2i4,4f18.10)') mg6,mm,mm,vector2(mg6),vector1(mm),wl,vector1(mm)

    !if (mg6 == ndr) then
    !  write(nf2,'(a9,3i6,3f20.14)') ' dbl_ext ',mg6,mm,mm,wl,vector1(mm),vector2(mg6)
    !  write(u6,'(a8,3i6,2f20.14)') ' dbl_ext ',mg1,mg2,mg3,wl,vector1(mm)
    !end if
end select

return

end subroutine prodel_2
