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

! ci diagonal elements
subroutine diagonal_loop_wyb()  !  for norb_act<>0

#include "drt_h.fh"
#include "intsort_h.fh"
#include "pl_structure_h.fh"
!common /casrst/ja(max_node),jb(max_node),jm(0:max_node),jj(4,0:max_node),kk(0:max_node),no(0:max_innorb),jv,jd(8),jt(8),js(8)
!common /sub_drt/jpad,jpae,ipae,ndim,nohy,ihy(max_wei),jj_sub(4,0:max_node),iy(4,0:max_node),jphy(max_node)

do lr0=2,norb_all
  do lr=1,lr0-1
    vdint(lr,lr0) = voint(lr0,lr)-vdint(lr0,lr)-vdint(lr0,lr)   ! 520
    !write(6,'(2i4,3f14.8)') lr,lr0,voint(lr0,lr),vdint(lr0,lr),vdint(lr,lr0)
  end do
end do
!write(6,*) '               ***** start h-diaelm *****'
vector1(1:nci_dim) = vpotnuc
!wl8 = hnil*(hnil-1)*vmd(lr,lr)*0.5d0+hnil*vo(lr,lr)   ! 800

ndimsum = 0
jpae = jv
ipae = 1
jaedownwei = iseg_downwei(ipae)
do jpad_=1,mxnode
  jpad = jpad_ ! jpad is in common block, is this necessary?
  iw_sta(jpad,ipae) = ndimsum
  if (nu_ad(jpad) == 0) cycle
  call seg_drt()
  iwupwei = jpad_upwei(jpad)
  iw_downwei(jpad,ipae) = ndim
  ndimsum = ndimsum+ndim*jaedownwei*iwupwei
  if (ndim == 0) cycle
  call diagonal_act_d()
  call diagonal_act_c()
end do
do im=1,ng_sm
  jpae = jd(im)
  ipae = 1+im
  if (nu_ae(ipae) == 0) cycle
  jaedownwei = iseg_downwei(ipae)
  do jpad_=1,mxnode
    jpad = jpad_ ! jpad is in common block, is this necessary?
    iw_sta(jpad,ipae) = ndimsum
    if (nu_ad(jpad) == 0) cycle
    call seg_drt()
    iwupwei = jpad_upwei(jpad)
    iw_downwei(jpad,ipae) = ndim
    !if (jpad >= 26) write(6,*)
    ndimsum = ndimsum+ndim*jaedownwei*iwupwei
    if (ndim == 0) cycle
    call diagonal_act_d()
    call diagonal_act_c()
  end do
end do
do im=1,ng_sm
  jpae = jt(im)
  ipae = 9+im
  if (nu_ae(ipae) == 0) cycle
  jaedownwei = iseg_downwei(ipae)
  do jpad_=1,mxnode
    jpad = jpad_ ! jpad is in common block, is this necessary?
    iw_sta(jpad,ipae) = ndimsum
    if (nu_ad(jpad) == 0) cycle
    call seg_drt()
    iwupwei = jpad_upwei(jpad)
    iw_downwei(jpad,ipae) = ndim
    ndimsum = ndimsum+ndim*jaedownwei*iwupwei
    if (ndim == 0) cycle
    call diagonal_act_d()
    call diagonal_act_c()
  end do
end do
do im=1,ng_sm
  jpae = js(im)
  ipae = 17+im
  jaedownwei = iseg_downwei(ipae)
  if (nu_ae(ipae) == 0) cycle
  do jpad_=1,mxnode
    jpad = jpad_ ! jpad is in common block, is this necessary?
    iw_sta(jpad,ipae) = ndimsum
    if (nu_ad(jpad) == 0) cycle
    call seg_drt()
    iwupwei = jpad_upwei(jpad)
    iw_downwei(jpad,ipae) = ndim
    ndimsum = ndimsum+ndim*jaedownwei*iwupwei
    if (ndim == 0) cycle
    call diagonal_act_d()
    call diagonal_act_c()
  end do
end do
call diagonal_dbl()
call diagonal_ext()
outer: do ipae_=1,25
  ipae = ipae_ ! ipae is in common block, is this necessary?
  do jpad_=1,mxnode
    jpad = jpad_ ! jpad is in common block, is this necessary?
    if (iw_sta(jpad,ipae) /= 0) then
      nci_h0 = iw_sta(jpad,ipae)
      exit outer
    end if
  end do
end do outer

!do i=1,nci_dim
!  write(6,'(i8,f18.8)') i,vector1(i)
!end do
!call abend()

return

end subroutine diagonal_loop_wyb

subroutine diagonal_act_c()

#include "drt_h.fh"
#include "intsort_h.fh"
#include "pl_structure_h.fh"
!common /casrst/ja(max_node),jb(max_node),jm(0:max_node),jj(4,0:max_node),kk(0:max_node),no(0:max_innorb),jv,jd(8),jt(8),js(8)
!common /sub_drt/jpad,jpae,ipae,ndim,nohy,ihy(max_wei),jj_sub(4,0:max_node),iy(4,0:max_node),jphy(max_node)
dimension ndr(max_innorb)
integer, pointer :: jph(:), jeh(:), jwh(:)
real*8, pointer :: th(:), thh(:)
#include "ptlph.fh"
#include "ptlphv.fh"
real*8, allocatable :: te(:), tee(:)
integer, allocatable :: jpe(:), jee(:), jwe(:)

!write(6,*) '               ***** start h-diaelm *****'
!write(6,*) jpad,jpae
allocate(te(maxpl),tee(maxpl),jpe(maxpl),jee(maxpl),jwe(maxpl))
ndr = 0
if (norb_act == 0) then
  mh = 1
  th(1) = 1.d0
  thh(1) = 1.d0
  call diagonal_link_dae(mh)
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
  if ((vlop0 == 0.0d0) .and. (vlop1 == 0.0d0)) cycle
  mpe = jj_sub(idl,jp)
  iwa = iy(idl,jp)
  call diagonal_link_ad(mpe,iwa,vlop0,vlop1)
end do
!***********************************************************************
!write(6,*) ad(i)
!***********************************************************************
lr = norb_dz+1

40 continue
if (lr == norb_inn) then
  if (mh /= 0) call diagonal_link_dae(mh)
  return
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
    !if (lr == 1) goto 20
    ! link loop
    call smidc2(isq,w,ww,mw,ind1,jeb)
    me = me+1
    jwe(me) = jwh(m)
    if (idl /= 1) jwe(me) = jwe(me)+iy(idl,je)
    jee(me) = jj_sub(idl,je)
    te(me) = th(m)*w
    tee(me) = thh(m)*ww
    jpe(me) = jph(m)
    ! complete a loop 'v'
    !20 continue
    if (ind1 == 1) cycle
    call stml(isq,w,ww,mw,ind1-1,jeb)
    vlop0 = th(m)*w
    vlop1 = thh(m)*ww
    if ((vlop0 == 0.0d0) .and. (vlop1 == 0.0d0)) cycle
    mp = mp+1
    mpe = jj_sub(idl,je)
    iwa = jwh(m)
    if (idl /= 1) iwa = iwa+iy(idl,je)
    call diagonal_link_ad(mpe,iwa,vlop0,vlop1)
    !*****   520  ******************************************************
  end do
end do
!***********************************************************************
!write(6,*) ad(i)
!***********************************************************************
do m=1,me
  th(m) = te(m)
  te(m) = 0.0d0
  thh(m) = tee(m)
  tee(m) = 0.0d0
  jwh(m) = jwe(m)
  jwe(m) = 0
  jeh(m) = jee(m)
  jee(m) = 0
  jph(m) = jpe(m)
  jpe(m) = 0
end do
mh = me
if (ndr(lr) < mh) ndr(lr) = mh
goto 40
do m=1,mh
  th(m) = 0.0d0
  thh(m) = 0.0d0
  jwh(m) = 0
  jph(m) = 0
  jeh(m) = 0
end do

return

end subroutine diagonal_act_c

subroutine diagonal_act_d()

#include "drt_h.fh"
#include "intsort_h.fh"
#include "pl_structure_h.fh"
!common /casrst/ja(max_node),jb(max_node),jm(0:max_node),jj(4,0:max_node),kk(0:max_node),no(0:max_innorb),jv,jd(8),jt(8),js(8)
!common /sub_drt/jpad,jpae,ipae,ndim,nohy,ihy(max_wei),jj_sub(4,0:max_node),iy(4,0:max_node),jphy(max_node)
dimension ndr(max_innorb)
integer, pointer :: jph(:), jeh(:), jwh(:)
real*8, pointer :: th(:), thh(:)
#include "ptlph.fh"
#include "ptlphv.fh"
real*8, allocatable :: te(:), tee(:)
integer, allocatable :: jpe(:), jee(:), jwe(:)

!write(6,*) '               ***** start h-diaelm *****'
!write(6,*) '   diagonal_act_d:',jpad,ipae
allocate(te(maxpl),tee(maxpl),jpe(maxpl),jee(maxpl),jwe(maxpl))
ndr = 0
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
      wt = voint(lr0,lr0)    ! hnil=1
      jw = iy(idl,jp)
      call prodel(3,wt,jp,mpe,jw)
    end do
    mpe = jj_sub(4,jp)
    if (mpe /= 0) then
      wt = vdint(lr0,lr0)+2.d0*voint(lr0,lr0)     !idl=4 hnil=2
      jw = iy(4,jp)
      call prodel(3,wt,jp,mpe,jw)
    end if
  end do
end do
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
  !write(6,*) ad(i)
  !*********************************************************************
  lr = lr0
  if (ndr(lr) < mh) ndr(lr) = mh
40 continue
  if (lr == norb_inn) then
    call diagonal_link_ae(mh)
    goto 5
  end if
  lr = lr+1
  me = 0
  mp = 0
  do m=1,mh
    je = jeh(m)
    jeb = jb(je)
    jp = jph(m)
    do idl=1,4
      if (jj_sub(idl,je) == 0) goto 17
      ind1 = idl
      if (lr == 1) goto 20
      ! link loop
      call smidc2(isq,w,ww,mw,ind1,jeb)
      me = me+1
      jwe(me) = jwh(m)
      if (idl /= 1) jwe(me) = jwe(me)+iy(idl,je)
      jee(me) = jj_sub(idl,je)
      te(me) = th(m)*w
      tee(me) = thh(m)*ww
      jpe(me) = jph(m)
      ! complete a loop 'v'
20    continue
      if (ind1 == 1) goto 17
      call stml(isq,w,ww,mw,ind1-1,jeb)
      vlop0 = th(m)*w
      vlop1 = thh(m)*ww
      if ((vlop0 == 0.0d0) .and. (vlop1 == 0.0d0)) goto 17
      mp = mp+1
      mpe = jj_sub(idl,je)
      iwa = jwh(m)
      if (idl /= 1) iwa = iy(idl,je)+iwa
      wt = (vlop0-vlop1)*voint(lr,lr0)-2.d0*vlop0*vdint(lr,lr0)

      call prodel(3,wt,jp,mpe,iwa)
      !*****   520  ****************************************************
17    continue
    end do
  end do
  !*********************************************************************
  !write(6,*) ad(i)
  !*********************************************************************
  do m=1,me
    th(m) = te(m)
    te(m) = 0.0d0
    thh(m) = tee(m)
    tee(m) = 0.0d0
    jwh(m) = jwe(m)
    jwe(m) = 0
    jeh(m) = jee(m)
    jee(m) = 0
    jph(m) = jpe(m)
    jpe(m) = 0
  end do
  mh = me
  if (ndr(lr) < mh) ndr(lr) = mh
  goto 40
  do m=1,mh
    th(m) = 0.0d0
    thh(m) = 0.0d0
    jwh(m) = 0
    jph(m) = 0
    jeh(m) = 0
  end do
5 continue
end do

return

end subroutine diagonal_act_d

subroutine diagonal_link_ae(mh)

#include "drt_h.fh"
#include "intsort_h.fh"
#include "pl_structure_h.fh"
!common /casrst/ja(max_node),jb(max_node),jm(0:max_node),jj(4,0:max_node),kk(0:max_node),no(0:max_innorb),jv,jd(8),jt(8),js(8)
!common /sub_drt/jpad,jpae,ipae,ndim,nohy,ihy(max_wei),jj_sub(4,0:max_node),iy(4,0:max_node),jphy(max_node)
integer, pointer :: jph(:), jeh(:), jwh(:)
real*8, pointer :: th(:), thh(:)
#include "ptlph.fh"
#include "ptlphv.fh"
data dsq2,vsq2/1.414213562373095d0,0.7071067811865d0/

!dsq2 = sqrt(2.d0)
!vsq2 = 1/sqrt(2.d0)
!dsq3vsq2 = sqrt(3.d0)/sqrt(2.d0)
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
  !write(6,*) 'ip,jpe,ind0',ip,jpe,ind0
  iwa = jwh(ip)
  vlop0 = th(ip)
  vlop1 = thh(ip)

  !wl5 = (vlop0-vlop1)*vo(lr0,lr)-2.0d0*vlop0*vmd(lr0,lr)
  !wl8 = vlop0*(vo(lr0,lr0)+(vlop0-1)*0.5*vmd(lr0,lr0))
  ! two-index,one-loop 520
  ! 0=<a,j,k,a>:13,14(ss=3),38(tt=2),50(dd=1)
  goto(100,200,300),ityae
  !link arc_d
100 continue
  !zz = '  g50  '
  wg50 = vlop0*vsq2
  wwg50 = -vlop1*sqrt(3.d0)/sqrt(2.d0)
  do la=1,norb_ext
    ma = lsm(la)
    if (ma /= imae) cycle
    lra = norb_all-la+1
    iwe = iwe+1
    wld = (wg50-wwg50)*voint(lra,lr0)-2.d0*wg50*vdint(lra,lr0)
    !write(6,'(a11,2i3,i6,1x,5f10.4)') zz,lr0,la,jwl,vo(lr0,la),vmd(lr0,la),wg50,wwg50,wl
    call prodel(5,wld,jp,iwa,iwe)
  end do
  goto 108
200 continue
  !zz = '  g38,39  '
  wg38 = -vlop0*vsq2
  wwg38 = vlop1
  do ima=1,ng_sm
    imb = mul_tab(ima,imae)
    if (imb > ima) cycle
    do la=ibsm_ext(ima),iesm_ext(ima)
      lra = norb_all-la+1
      lbsta = ibsm_ext(imb)
      lbend = iesm_ext(imb)
      if (ima == imb) lbend = la-1
      do lb=lbsta,lbend
        lrb = norb_all-lb+1
        iwe = iwe+1
        wlt = (wg38-wwg38)*(voint(lra,lr0)+voint(lrb,lr0))-2.d0*wg38*(vdint(lra,lr0)+vdint(lrb,lr0))
        !write(6,*) ' 520 r0,la,lb ',vo(r0,la),vo(r0,lb),vmd(r0,la),vmd(r0,lb)

        call prodel(5,wlt,jp,iwa,iwe)
      end do
    end do
  end do
  goto 108
300 continue
  !zz = '  g14,15  '
  wg14 = -vlop0*vsq2
  do ima=1,ng_sm
    imb = mul_tab(ima,imae)
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
        wls = wg14*(voint(lra,lr0)+voint(lrb,lr0))-2.d0*wg14*(vdint(lra,lr0)+vdint(lrb,lr0))
        call prodel(5,wls,jp,iwa,iwe)
      end do
    end do
  end do
  if (ipae /= 18) goto 108
  !zz = '  g13     '
  wg13 = -vlop0*dsq2
  do la=1,norb_ext
    lra = norb_all-la+1
    iwe = iwe+1
    wls = wg13*(voint(lra,lr0)-2.d0*vdint(lra,lr0))
    call prodel(5,wls,jp,iwa,iwe)
    !write(6,*) ' g13 ',vo(lr0,la),vo(lr0,lb),vmd(lr0,la),vmd(lr0,lb)
  end do
108 continue
end do
mh = 0

return

end subroutine diagonal_link_ae

subroutine diagonal_link_ad(mpe,iwa,vlop0,vlop1)

#include "drt_h.fh"
#include "intsort_h.fh"
#include "pl_structure_h.fh"
!common /casrst/ja(max_node),jb(max_node),jm(0:max_node),jj(4,0:max_node),kk(0:max_node),no(0:max_innorb),jv,jd(8),jt(8),js(8)
!common /sub_drt/jpad,jpae,ipae,ndim,nohy,ihy(max_wei),jj_sub(4,0:max_node),iy(4,0:max_node),jphy(max_node)
integer, pointer :: jph(:), jeh(:), jwh(:)
real*8, pointer :: th(:), thh(:)
#include "ptlph.fh"
#include "ptlphv.fh"
data dsq2,vsq2/1.414213562373095d0,0.7071067811865d0/

!dsq2 = sqrt(2.d0)
!vsq2 = 1/sqrt(2.d0)
!dsq3vsq2 = sqrt(3.d0)/sqrt(2.d0)
if (norb_dz == 0) return
ityad = 1
if (jpad /= 1) ityad = (jpad-1)/8+2
imad = mod(jpad-1,8)
if (imad == 0) then
  ityad = ityad-1
  imad = 8
end if
lra = kk(mpe)-1
goto(100,200,300,400,500,600),ityad
! v: d&r&l(3)
100 continue
fqi = -fg
vl0 = fqi*dsq2*vlop0
wlv = 0.d0
!do lr=norb_frz+1,norb_dz
do lr=1,norb_dz
  wlv = wlv-vl0*vdint(lr,lra)
end do
call prodel(4,wlv,mpe,0,iwa)
return
!jpad = jd(im)
200 continue
fqi = -fg
do lri=norb_frz+1,norb_dz
  imd = mul_tab(lsm_inn(lri),ns_sm)
  if (imd /= imad) cycle
  iwd = jud(lri)

  ! d: d&r&l(2)
  vl0 = fqi*vsq2*vlop0
  vl1 = pd*vlop1
  wld = -2.0d0*vl0*vdint(lra,lri)+(vl0-vl1)*voint(lra,lri)
  ! d: d&r&l(3)+c"(2)
  vl0 = fqi*dsq2*vlop0
  do lr=1,lri-1
    wld = wld+vl0*(voint(lra,lr)-2.d0*vdint(lra,lr))
  end do
  ! d: d&r&l(3)
  do lr=lri+1,norb_dz
    wld = wld+vl0*(voint(lra,lr)-2.d0*vdint(lra,lr))
  end do

  call prodel(4,wld,mpe,iwd,iwa)
end do
return
!jpad = jt(im)
300 continue
fqi = fg
iwt = 0
do lri=norb_frz+1,norb_dz
  imi = mul_tab(lsm_inn(lri),ns_sm)
  do lrj=lri+1,norb_dz
    imj = lsm_inn(lrj)
    imij = mul_tab(imi,imj)
    if (imij /= imad) cycle
    iwt = just(lri,lrj)

    ! t: d&r&l(2)
    ! t: d&r&l(2)+c"(2)
    vl0 = fqi*vsq2*vlop0
    vl1 = pt*vlop1
    wlt = -2.0d0*vl0*vdint(lra,lri)+(vl0-vl1)*voint(lra,lri)
    wlt = wlt-2.0d0*vl0*vdint(lra,lrj)+(vl0-vl1)*voint(lra,lrj)
    vl0 = fqi*dsq2*vlop0
    ! t: d&r&l(3)+c"(2)+c"(2)
    do lr=1,lri-1
      wlt = wlt+vl0*vdint(lr,lra)
    end do
    ! t: d&r&l(3)+c"(2)
    do lr=lri+1,lrj-1
      wlt = wlt+vl0*vdint(lr,lra)
    end do
    ! t: d&r&l(3)
    do lr=lrj+1,norb_dz
      wlt = wlt+vl0*vdint(lr,lra)
    end do
    call prodel(4,wlt,mpe,iwt,iwa)
  end do
end do
return
400 continue
fqi = fg
iws = 0
do lri=norb_frz+1,norb_dz
  if (imad == ns_sm) then
    lrj = lri
    iws = just(lri,lri)
    ! s: d&r&l(3)+c"(0)
    ! s: d&r&l(3)
    vl0 = fqi*dsq2*vlop0
    wls = 0.d0
    do lr=1,norb_dz
      if (lr == lri) cycle
      wls = wls+vl0*(voint(lra,lr)-2.d0*vdint(lra,lr))
    end do
    call prodel(4,wls,mpe,iws,iwa)
  end if

  imi = mul_tab(lsm_inn(lri),ns_sm)
  lrjsta = lri+1
  do lrj=lrjsta,norb_dz
    imj = lsm_inn(lrj)
    imij = mul_tab(imi,imj)
    if (imij /= imad) cycle
    iws = just(lri,lrj)
    ! s1: d&r&l(1)
    vl0 = fqi*vsq2*vlop0
    vl1 = ps1*vlop1
    wls = -2.0d0*vl0*vdint(lra,lrj)+(vl0-vl1)*voint(lra,lrj)
    ! s4: d&r&l(2)+c"(1)
    vl0 = fqi*vsq2*vlop0
    vl1 = ps4*vlop1
    wls = wls-2.0d0*vl0*vdint(lra,lri)+(vl0-vl1)*voint(lra,lri)

    vl0 = fqi*dsq2*vlop0
    ! s: d&r&l(3)+c"(2)+c"(1)
    do lr=1,lri-1
      wls = wls+vl0*vdint(lr,lra)
    end do
    ! s: d&r&l(3)+c"(1)
    do lr=lri+1,lrj-1
      wls = wls+vl0*vdint(lr,lra)
    end do
    ! s: d&r&l(3)
    do lr=lrj+1,norb_dz
      wls = wls+vl0*vdint(lr,lra)
    end do
    call prodel(4,wls,mpe,iws,iwa)
  end do
end do

if (jb_sys == 0) return      !any difference when jb_sys=1 and jb_sy
fqi = fg
do lri=norb_frz+1,norb_dz
  imi = mul_tab(lsm_inn(lri),ns_sm)
  do lrj=lri+1,norb_dz
    imj = lsm_inn(lrj)
    imij = mul_tab(imi,imj)
    if (imij /= imad) cycle
    !iws = iws+1
    iws = just(lrj,lri)
    ! s1: d&r&l(1)-c"(2)
    vl0 = fqi*vsq2*vlop0
    vl1 = ps3*vlop1
    wls = -2.0d0*vl0*vdint(lra,lri)+(vl0-vl1)*voint(lra,lri)

    ! s3: (11)d&r&l(2)
    vl0 = fqi*vsq2*vlop0
    vl1 = ps2*vlop1
    wls = wls-2.0d0*vl0*vdint(lra,lrj)+(vl0-vl1)*voint(lra,lrj)

    vl0 = fqi*dsq2*vlop0
    ! s: d&r&l(3)+c"(1)+c"(2)
    do lr=1,lri-1
      wls = wls+vl0*vdint(lr,lra)
    end do
    ! s: d&r&l(3)+c"(2)
    do lr=lri+1,lrj-1
      wls = wls+vl0*vdint(lr,lra)
    end do
    ! s: d&r&l(3)
    do lr=lrj+1,norb_dz
      wls = wls+vl0*vdint(lr,lra)
    end do
    call prodel(4,wls,mpe,iws,iwa)
  end do
end do
return

500 continue
fqi = -fg
iwd = 0

do lri=norb_frz+1,norb_dz
  imd = mul_tab(lsm_inn(lri),ns_sm)
  if (imd /= imad) cycle
  iwd = jud(lri)

  ! dd1: d&r&l(1)
  vl0 = fqi*vsq2*vlop0
  vl1 = pdd*vlop1
  !if ((mpe == 43) .and. (iwd == 1.)) then
  !  write(6,*)
  !  write(6,*) lri,lra,vl0,vl1,wld,vlop0,vlop1
  !end if

  wld = -2.0d0*vl0*vdint(lra,lri)+(vl0-vl1)*voint(lra,lri)

  vl0 = fqi*dsq2*vlop0

  ! d: d&r&l(3)+c"(1)
  do lr=1,lri-1
    wld = wld+vl0*vdint(lr,lra)
  end do

  ! d: d&r&l(3)
  do lr=lri+1,norb_dz
    wld = wld+vl0*vdint(lr,lra)
  end do

  call prodel(4,wld,mpe,iwd,iwa)
end do
return

600 continue
fqi = fg
iwt = 0
do lri=norb_frz+1,norb_dz
  imi = mul_tab(lsm_inn(lri),ns_sm)
  do lrj=lri+1,norb_dz
    imj = lsm_inn(lrj)
    imij = mul_tab(imi,imj)
    if (imij /= imad) cycle
    iwt = just(lri,lrj)

    ! tt: d&r&l(1)
    ! tt: d&r&l(1)+c"(1)
    vl0 = fqi*vsq2*vlop0
    vl1 = ptt*vlop1
    wlt = -2.0d0*vl0*vdint(lra,lri)+(vl0-vl1)*voint(lra,lri)
    wlt = wlt-2.0d0*vl0*vdint(lra,lrj)+(vl0-vl1)*voint(lra,lrj)
    vl0 = fqi*dsq2*vlop0
    ! t: d&r&l(3)+c"(1)+c"(1)
    do lr=1,lri-1
      wlt = wlt+vl0*vdint(lr,lra)
    end do
    ! t: d&r&l(3)+c"(1)
    do lr=lri+1,lrj-1
      wlt = wlt+vl0*vdint(lr,lra)
    end do
    ! t: d&r&l(3)
    do lr=lrj+1,norb_dz
      wlt = wlt+vl0*vdint(lr,lra)
    end do

    call prodel(4,wlt,mpe,iwt,iwa)
  end do
end do

return

end subroutine diagonal_link_ad

subroutine diagonal_link_dae(mh)

#include "drt_h.fh"
#include "intsort_h.fh"
#include "pl_structure_h.fh"
!common /casrst/ja(max_node),jb(max_node),jm(0:max_node),jj(4,0:max_node),kk(0:max_node),no(0:max_innorb),jv,jd(8),jt(8),js(8)
!common /sub_drt/jpad,jpae,ipae,ndim,nohy,ihy(max_wei),jj_sub(4,0:max_node),iy(4,0:max_node),jphy(max_node)
integer, pointer :: jph(:), jeh(:), jwh(:)
real*8, pointer :: th(:), thh(:)
#include "ptlph.fh"
#include "ptlphv.fh"
data dsq2,vsq2/1.414213562373095d0,0.7071067811865d0/

!dsq2 = sqrt(2.d0)
!vsq2 = 1/sqrt(2.d0)
!dsq3vsq2 = sqrt(3.d0)/sqrt(2.d0)
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

  goto(100,200,300,400,500,600),ityad

  !jpad = 1
100 continue
  if (abs(vlop0) < 1e-30) goto 108
  fqi = fg
  ! v: d&r&l(3)
  lri = 0
  lrj = 0
  iwd = 0
  vij0 = 0.d0
  vij1 = 0.d0
  vij2 = 0.d0
  vl0 = fqi*dsq2*vlop0
  call diagonal_call_dae(lri,lrj,iwd,iwa,vij0,vij1,vij2,vl0)
  goto 108
  !jpad = jd(im)
200 continue
  fqi = -fg
  lrj = 0
  do lri=norb_frz+1,norb_dz
    imd = mul_tab(lsm_inn(lri),ns_sm)
    if (imd /= imad) cycle
    iwd = jud(lri)

    ! d: d&r&l(2)
    vij0 = fqi*vsq2*vlop0
    vij1 = pd*vlop1
    vij2 = 0.d0
    ! d: d&r&l(3)+c"(2)
    ! d: d&r&l(3)
    vl0 = fqi*dsq2*vlop0

    call diagonal_call_dae(lri,lrj,iwd,iwa,vij0,vij1,vij2,vl0)

  end do
  goto 108
  !jpad = jt(im)
300 continue
  fqi = fg
  iwt = 0
  do lri=norb_frz+1,norb_dz
    imi = mul_tab(lsm_inn(lri),ns_sm)
    do lrj=lri+1,norb_dz
      imj = lsm_inn(lrj)
      imij = mul_tab(imi,imj)
      if (imij /= imad) cycle
      iwt = just(lri,lrj)

      ! t: d&r&l(2)
      ! t: d&r&l(2)+c"(2)
      vij0 = fqi*vsq2*vlop0
      vij1 = pt*vlop1
      vij2 = vij1
      ! t: d&r&l(3)+c"(2)+c"(2)
      ! t: d&r&l(3)+c"(2)
      ! t: d&r&l(3)
      vl0 = fqi*dsq2*vlop0

      call diagonal_call_dae(lri,lrj,iwt,iwa,vij0,vij1,vij2,vl0)
    end do
  end do
  goto 108
400 continue
  fqi = fg
  iws = 0
  do lri=norb_frz+1,norb_dz           !cc
    if (imad == ns_sm) then
      lrj = lri
      iws = just(lri,lri)
      vij0 = 0.d0
      vij1 = 0.d0
      vij2 = 0.d0
      ! s: d&r&l(3)+c"(0)
      ! s: d&r&l(3)
      vl0 = fqi*dsq2*vlop0
      call diagonal_call_dae(lri,lrj,iws,iwa,vij0,vij1,vij2,vl0)
    end if
    imi = mul_tab(lsm_inn(lri),ns_sm)
    lrjsta = lri+1
    do lrj=lrjsta,norb_dz
      imj = lsm_inn(lrj)
      imij = mul_tab(imi,imj)
      if (imij /= imad) cycle
      iws = just(lri,lrj)
      ! s2: d&r&l(2)
      ! s4: d&r&l(2)+c"(1)
      vij0 = fqi*vsq2*vlop0
      vij1 = ps1*vlop1
      vij2 = ps4*vlop1
      ! s: d&r&l(3)+c"(2)+c"(1)
      ! s: d&r&l(3)+c"(1)
      ! s: d&r&l(3)
      vl0 = fqi*dsq2*vlop0
      call diagonal_call_dae(lri,lrj,iws,iwa,vij0,vij1,vij2,vl0)
    end do
  end do
  if (jb_sys == 0) goto 108
  fqi = fg
  do lri=norb_frz+1,norb_dz
    imi = mul_tab(lsm_inn(lri),ns_sm)
    if (imad == ns_sm) lrjsta = lri
    do lrj=lri+1,norb_dz
      imj = lsm_inn(lrj)
      imij = mul_tab(imi,imj)
      if (imij /= imad) cycle
      iws = just(lrj,lri)
      ! s1: d&r&l(1)
      ! s3: d&r&l(1)+c"(2)
      vij0 = fqi*vsq2*vlop0
      vij1 = ps2*vlop1
      vij2 = ps3*vlop1
      ! s: d&r&l(3)+c"(1)+c"(2)
      ! s: d&r&l(3)+c"(2)
      ! s: d&r&l(3)
      vl0 = fqi*dsq2*vlop0
      call diagonal_call_dae(lri,lrj,iws,iwa,vij0,vij1,vij2,vl0)
    end do
  end do
  goto 108
500 continue
  fqi = -fg
  lrj = 0
  do lri=norb_frz+1,norb_dz
    imd = mul_tab(lsm_inn(lri),ns_sm)
    if (imd /= imad) cycle
    iwd = jud(lri)

    ! dd1: d&r&l(1)
    vij0 = fqi*vsq2*vlop0
    vij1 = pdd*vlop1
    vij2 = 0.d0
    ! d: d&r&l(3)+c"(1)
    ! d: d&r&l(3)
    vl0 = fqi*dsq2*vlop0

    call diagonal_call_dae(lri,lrj,iwd,iwa,vij0,vij1,vij2,vl0)

  end do
  goto 108

600 continue
  fqi = fg               !aa
  iwt = 0
  do lri=norb_frz+1,norb_dz
    imi = mul_tab(lsm_inn(lri),ns_sm)
    do lrj=lri+1,norb_dz
      imj = lsm_inn(lrj)
      imij = mul_tab(imi,imj)
      if (imij /= imad) cycle
      iwt = just(lri,lrj)

      vij0 = fqi*vsq2*vlop0
      ! tt: d&r&l(1)
      vij1 = ptt*vlop1
      ! tt: d&r&l(1)+c"(1)
      vij2 = ptt*vlop1
      ! t: d&r&l(3)+c"(1)+c"(1)
      ! t: d&r&l(3)+c"(1)
      ! t: d&r&l(3)
      vl0 = fqi*dsq2*vlop0
      call diagonal_call_dae(lri,lrj,iwt,iwa,vij0,vij1,vij2,vl0)
    end do
  end do

108 continue
end do

return

end subroutine diagonal_link_dae

subroutine diagonal_call_dae(lri,lrj,iwd,iwa,vij0,vij1,vij2,vl0)

#include "drt_h.fh"
#include "intsort_h.fh"
#include "pl_structure_h.fh"
!common /casrst/ja(max_node),jb(max_node),jm(0:max_node),jj(4,0:max_node),kk(0:max_node),no(0:max_innorb),jv,jd(8),jt(8),js(8)
!common /sub_drt/jpad,jpae,ipae,ndim,nohy,ihy(max_wei),jj_sub(4,0:max_node),iy(4,0:max_node),jphy(max_node)
integer, pointer :: jph(:), jeh(:), jwh(:)
real*8, pointer :: th(:), thh(:)
#include "ptlph.fh"
#include "ptlphv.fh"
data dsq2,vsq2/1.414213562373095d0,0.7071067811865d0/
data dsq3vsq2/1.224744871392d0/

!dsq2 = sqrt(2.d0)
!vsq2 = 1/sqrt(2.d0)
!dsq3vsq2 = sqrt(3.d0)/sqrt(2.d0)
if (norb_dz == 0) return
if (ipae == 1) return   !could not link
ityae = (ipae-1)/8+1
imae = mod(ipae-1,8)
if (imae == 0) then
  ityae = ityae-1
  imae = 8
end if
iwe = 0
! 520=<a,j,k,a>:13,14(ss=3),38(tt=2),50(dd=1)
goto(100,200,300),ityae
! link arc_d
100 continue
!zz = '  g50  '
do la=ibsm_ext(imae),iesm_ext(imae)
  lra = norb_all-la+1
  iwe = iwe+1
  vlop0 = vl0*vsq2
  wl = 0.d0
  do lr=1,norb_dz
    if (lr == lri) cycle
    if (lr == lrj) cycle
    wl = wl+vlop0*vdint(lr,lra)    !db space drl(33)- ext space -
  end do
  if (lri >= lrj) then
    vlop0 = vij0*vsq2
    vlop1 = -vij1*dsq3vsq2
    wl = wl+(vlop0-vlop1)*voint(lra,lri)-2.d0*vlop0*vdint(lra,lri)
    vlop1 = -vij2*dsq3vsq2
    wl = wl+(vlop0-vlop1)*voint(lra,lrj)-2.d0*vlop0*vdint(lra,lrj)
    !write(6,'(a11,2i3,i6,1x,5f10.4)') zz,lr0,la,jwl,vo(lr0,la),vmd(lr0,la),wg50,wwg50,wl
  else
    vlop0 = vij0*vsq2
    vlop1 = -vij1*dsq3vsq2         !db space (22)drl(11)- ext space -g
    wl = wl+(vlop0-vlop1)*voint(lra,lrj)-2.d0*vlop0*vdint(lra,lrj)
    vlop1 = -vij2*dsq3vsq2         !db space drl(22)c"(11)- ext space
    wl = wl+(vlop0-vlop1)*voint(lra,lri)-2.d0*vlop0*vdint(lra,lri)
  end if
  call prodel(6,wl,iwd,iwa,iwe)
end do
goto 108
200 continue
!zz = '  g38,39  '
do ima=1,ng_sm
  imb = mul_tab(ima,imae)
  if (imb > ima) cycle
  do la=ibsm_ext(ima),iesm_ext(ima)
    lra = norb_all-la+1
    lbsta = ibsm_ext(imb)
    lbend = iesm_ext(imb)
    if (ima == imb) lbend = la-1
    do lb=lbsta,lbend
      lrb = norb_all-lb+1
      iwe = iwe+1
      volalb = 0.d0
      vd2lalb = 0.d0
      do lr=1,norb_dz
        if (lr == lri) cycle
        if (lr == lrj) cycle
        volalb = volalb+(voint(lra,lr)+voint(lrb,lr))
        vd2lalb = vd2lalb-2.d0*(vdint(lra,lr)+vdint(lrb,lr))
      end do

      vlop0 = -vl0*vsq2
      wl = vlop0*(volalb+vd2lalb)
      if (lri >= lrj) then
        vlop0 = -vij0*vsq2
        vlop1 = vij1
        wl = wl+(vlop0-vlop1)*(voint(lra,lri)+voint(lrb,lri))-2.d0*vlop0*(vdint(lra,lri)+vdint(lrb,lri))
        vlop1 = vij2
        wl = wl+(vlop0-vlop1)*(voint(lra,lrj)+voint(lrb,lrj))-2.d0*vlop0*(vdint(lra,lrj)+vdint(lrb,lrj))
      else
        vlop0 = -vij0*vsq2
        vlop1 = vij1
        wl = wl+(vlop0-vlop1)*(voint(lra,lrj)+voint(lrb,lrj))-2.d0*vlop0*(vdint(lra,lrj)+vdint(lrb,lrj))
        vlop1 = vij2
        wl = wl+(vlop0-vlop1)*(voint(lra,lri)+voint(lrb,lri))-2.d0*vlop0*(vdint(lra,lri)+vdint(lrb,lri))
      end if
      !write(6,*) ' 520 r0,la,lb ',vo(r0,la),vo(r0,lb),vmd(r0,la),vmd(r0,lb)

      call prodel(6,wl,iwd,iwa,iwe)
    end do
  end do
end do
goto 108
300 continue
!zz = '  g14,15  '
do ima=1,ng_sm
  imb = mul_tab(ima,imae)
  if (imb > ima) cycle
  do la=ibsm_ext(ima),iesm_ext(ima)
    lra = norb_all-la+1
    lbsta = ibsm_ext(imb)
    lbend = iesm_ext(imb)
    if (ima == imb) lbend = la-1
    do lb=lbsta,lbend
      lrb = norb_all-lb+1
      iwe = iwe+1
      volalb = 0.d0
      vd2lalb = 0.d0
      do lr=1,norb_dz
        if (lr == lri) cycle
        if (lr == lrj) cycle
        volalb = volalb+(voint(lra,lr)+voint(lrb,lr))
        vd2lalb = vd2lalb-2.d0*(vdint(lra,lr)+vdint(lrb,lr))
      end do

      wg14 = -vl0*vsq2
      wl = wg14*(volalb+vd2lalb)

      if ((jpad == 18) .and. (lri == lrj)) goto 301
      wg14 = -vij0*vsq2
      wl = wl+wg14*(voint(lra,lri)+voint(lrb,lri))-2.d0*wg14*(vdint(lra,lri)+vdint(lrb,lri))
      wl = wl+wg14*(voint(lra,lrj)+voint(lrb,lrj))-2.d0*wg14*(vdint(lra,lrj)+vdint(lrb,lrj))
      !write(6,*) ' 520 r0,la,lb ',vo(r0,la),vo(r0,lb),vmd(r0,la),vmd(r0,lb)
301   continue
      call prodel(6,wl,iwd,iwa,iwe)
    end do
  end do
end do

if (ipae /= 18) return
!zz = '  g13     '
do la=1,norb_ext
  lra = norb_all-la+1
  iwe = iwe+1
  vovdla = 0.d0
  do lr=1,norb_dz
    if (lr == lri) cycle
    if (lr == lrj) cycle
    vovdla = vovdla+vdint(lr,lra)
  end do

  wg13 = -vl0*dsq2
  wl = wg13*vovdla
  wg13 = -vij0*dsq2
  wl = wl+wg13*(vdint(lri,lra)+vdint(lrj,lra))

  call prodel(6,wl,iwd,iwa,iwe)
  !write(6,*) ' g13 ',vo(lr0,la),vo(lr0,lb),vmd(lr0,la),vmd(lr0,lb)
end do

108 continue
return

end subroutine diagonal_call_dae

subroutine diagonal_dbl()

#include "drt_h.fh"
#include "intsort_h.fh"
#include "pl_structure_h.fh"
!common /casrst/ja(max_node),jb(max_node),jm(0:max_node),jj(4,0:max_node),kk(0:max_node),no(0:max_innorb),jv,jd(8),jt(8),js(8)
!common /sub_drt/jpad,jpae,ipae,ndim,nohy,ihy(max_wei),jj_sub(4,0:max_node),iy(4,0:max_node),jphy(max_node)
data dzero/0.d0/

if (norb_dz == 0) return
wt0 = dzero
do lr=1,norb_dz
  wt0 = wt0+voint(lr,lr)+voint(lr,lr)+vdint(lr,lr)
end do
do ipae_=1,25
  ipae = ipae_ ! ipae is in common block, is this necessary?
  if (nu_ae(ipae) == 0) cycle
  iwdownv = iw_downwei(1,ipae)
  do iwa=0,iwdownv-1
    !zz = ' doub_800_v'
    iwad = iwalk_ad(1,ipae,iwa,0)
    call prodel(1,wt0,0,ipae,iwad)
    !zz = ' doub_800_s'
  end do
end do
!jps = js(1)
do lr0=norb_frz+1,norb_dz
  mr0 = mul_tab(lsm_inn(lr0),ns_sm)
  iwd = jud(lr0)
  ! d_800
  jpad = 1+mr0
  jpad1 = jpad+24
  wld = wt0-voint(lr0,lr0)-vdint(lr0,lr0)
  wls = wld-voint(lr0,lr0)
  do ipae_=1,25
    ipae = ipae_ ! ipae is in common block, is this necessary?
    if (nu_ae(ipae) == 0) cycle
    iwdownv = iw_downwei(jpad,ipae)
    do iwa=0,iwdownv-1
      iwad = iwalk_ad(jpad,ipae,iwa,iwd)
      call prodel(1,wld,0,ipae,iwad)
    end do
  end do

  if (jb_sys > 0) then
    do ipae_=1,25
      ipae = ipae_ ! ipae is in common block, is this necessary?
      if (nu_ae(ipae) == 0) cycle
      iwdownv = iw_downwei(jpad1,ipae)
      do iwa=0,iwdownv-1
        iwad = iwalk_ad(jpad1,ipae,iwa,iwd)
        call prodel(1,wld,0,ipae,iwad)
      end do
    end do
  end if
  ! d_800
  jpad = 17+ns_sm
  iwd = just(lr0,lr0)
  do ipae_=1,25
    ipae = ipae_ ! ipae is in common block, is this necessary?
    if (nu_ae(ipae) == 0) cycle
    iwdownv = iw_downwei(jpad,ipae)
    do iwa=0,iwdownv-1
      iwad = iwalk_ad(jpad,ipae,iwa,iwd)
      call prodel(1,wls,0,ipae,iwad)
    end do
  end do

  if (lr0 == norb_dz) goto 100
  wld0 = wld

  do lr=lr0+1,norb_dz
    mr = mul_tab(mr0,lsm_inn(lr))
    jpat = 9+mr
    jpas = 17+mr
    jpat1 = jpat+24
    iws = just(lr0,lr)
    iwt = iws
    wld = wld0-voint(lr,lr)-vdint(lr,lr)
    do ipae_=1,25
      ipae = ipae_ ! ipae is in common block, is this necessary?
      if (nu_ae(ipae) == 0) cycle
      iwdownv = iw_downwei(jpat,ipae)
      do iwa=0,iwdownv-1
        iwad = iwalk_ad(jpat,ipae,iwa,iwt)
        call prodel(1,wld,0,ipae,iwad)
      end do
      if (jb_sys > 1) then
        iwdownv = iw_downwei(jpat1,ipae)
        do iwa=0,iwdownv-1
          iwad = iwalk_ad(jpat1,ipae,iwa,iwt)     !t1
          call prodel(1,wld,0,ipae,iwad)
        end do
      end if
      iwdownv = iw_downwei(jpas,ipae)
      do iwa=0,iwdownv-1
        iwad = iwalk_ad(jpas,ipae,iwa,iws)
        call prodel(1,wld,0,ipae,iwad)
      end do
      if (jb_sys > 0) then
        iws1 = just(lr,lr0)
        do iwa=0,iwdownv-1
          iwad = iwalk_ad(jpas,ipae,iwa,iws1)
          call prodel(1,wld,0,ipae,iwad)
        end do
      end if
    end do
  end do
100 continue
end do
!wl8 = 1/2*hnil*(hnil-1)*vmd(lr,lr)+hnil*voint(lr,lr)  !800
!wl5 = (vlop0-vlop1)*vo(lr0,lr)-2*vlop0*vmd(lr0,lr)
!wl5 = vlop0*vmd(lr,lr0)-vlop1*vo(lr0,lr)              !2000.11.26
wt0 = dzero
do lr0=2,norb_dz
  do lr=1,lr0-1
    ! vlop0=-2 vlop1=0
    wt0 = wt0-2.d0*vdint(lr,lr0)
  end do
end do
jpad = 1
iwd = 0
do ipae_=1,25
  ipae = ipae_ ! ipae is in common block, is this necessary?
  if (nu_ae(ipae) == 0) cycle
  iwdownv = iw_downwei(jpad,ipae)
  do iwa=0,iwdownv-1
    iwad = iwalk_ad(jpad,ipae,iwa,iwd)
    call prodel(1,wt0,0,ipae,iwad)
  end do
end do
do lrm=norb_frz+1,norb_dz
  mrm = mul_tab(lsm_inn(lrm),ns_sm)
  iws = just(lrm,lrm)
  iwd = jud(lrm)
  jpad = 1+mrm
  jpad1 = jpad+24
  jpas = 17+ns_sm
  ! d_520
  wld = wt0
  do lr=1,lrm-1
    wld = wld+vdint(lr,lrm)
  end do
  do lr0=lrm+1,norb_dz
    wld = wld+vdint(lrm,lr0)
  end do
  wls = wld
  do lr=lrm+1,norb_dz
    wls = wls+vdint(lrm,lr)
  end do
  do lr0=1,lrm-1
    wls = wls+vdint(lr0,lrm)
  end do
  do ipae_=1,25
    ipae = ipae_ ! ipae is in common block, is this necessary?
    if (nu_ae(ipae) == 0) cycle
    iwdownv = iw_downwei(jpad,ipae)
    do iwa=0,iwdownv-1
      iwad = iwalk_ad(jpad,ipae,iwa,iwd)
      call prodel(1,wld,0,ipae,iwad)
    end do
    if (jb_sys > 0) then
      iwdownv = iw_downwei(jpad1,ipae)
      do iwa=0,iwdownv-1
        iwad = iwalk_ad(jpad1,ipae,iwa,iwd)
        call prodel(1,wld,0,ipae,iwad)
      end do
    end if
    ! 520
    iwdownv = iw_downwei(jpas,ipae)
    do iwa=0,iwdownv-1
      iwad = iwalk_ad(jpas,ipae,iwa,iws)
      call prodel(1,wls,0,ipae,iwad)
    end do
  end do
end do
! 520
do lr0=norb_frz+1,norb_dz-1
  mr0 = mul_tab(lsm_inn(lr0),ns_sm)
  do lr=lr0+1,norb_dz
    mr = mul_tab(mr0,lsm_inn(lr))
    jpat = 9+mr
    jpas = 17+mr
    jpat1 = jpat+24
    iws = just(lr0,lr)
    iwt = iws
    iws1 = just(lr,lr0)
    if (jb_sys == 0) then
      wls = wt0+3.d0*(voint(lr,lr0)-vdint(lr,lr0))
      wlt = wt0+voint(lr,lr0)-3.d0*vdint(lr,lr0)
    end if
    ! 2*vo(lr0,lr)-4*vmd(lr0,lr)  w0=-2    w1=0
    !   vo(lr0,lr)+  vmd(lr0,lr)  w0=-1/2  w1=-3/2
    ! 2*vo(lr0,lr)-4*vmd(lr0,lr)  w0=-2    w1=0
    ! - vo(lr0,lr)+  vmd(lr0,lr)  w0=-1/2  w1=1/2
    if (jb_sys > 0) then
      db = jb_sys
      w1 = -(db+3)/(2.d0*db+2)
      wls = wt0+(1.50d0-w1)*voint(lr,lr0)-3.d0*vdint(lr,lr0)
      wlt = wt0+voint(lr,lr0)-3.d0*vdint(lr,lr0)
      w1 = -(db-1)/(2.d0*db+2)
      wls1 = wt0+(1.50d0-w1)*voint(lr,lr0)-3.d0*vdint(lr,lr0)
      ! 2*vo(lr0,lr)-4*vmd(lr0,lr)  w0=-2  w1=0
      ! (-1/2-w1)*vo(lr0,lr)+vmd(lr0,lr)
      ! w0=-1/2    w1=-(db-1)/(2*db+2)
    end if
    !if (jb_sys > 1) then
    !  wlt1 = wt0+voint(lr,lr0)-3*vdint(lr,lr0)
    !end if
    ! 2*vo(lr0,lr)-4*vmd(lr0,lr)  w0=-2  w1=0
    !  -vo(lr0,lr)+  vmd(lr0,lr)
    do lrg=1,lr0-1
      wls = wls+vdint(lrg,lr0)+vdint(lrg,lr)
      wlt = wlt+vdint(lrg,lr0)+vdint(lrg,lr)
      if (jb_sys > 0) then
        wls1 = wls1+vdint(lrg,lr0)+vdint(lrg,lr)
      end if
      !if (jb_sys > 1) then
      !  wlt1 = wlt1+vdint(lrg,lr0)+vdint(lrg,lr)
      !end if
      ! 2*vo(lr0,lr)-4*vmd(lr0,lr)  w0=-2  w1=0
      ! - vo(lr0,lr)+2*vmd(lr0,lr)  w0=-1  w1=0
    end do
    do lrg=lr0+1,lr-1
      wls = wls+vdint(lr0,lrg)+vdint(lrg,lr)
      wlt = wlt+vdint(lr0,lrg)+vdint(lrg,lr)
      if (jb_sys > 0) then
        wls1 = wls1+vdint(lr0,lrg)+vdint(lrg,lr)
      end if
      !if (jb_sys > 1) then
      !  wlt1 = wlt1+vdint(lr0,lrg)+vdint(lrg,lr)
      !end if
      ! 2*vo(lr0,lr)-4*vmd(lr0,lr)  w0=-2  w1=0
      ! - vo(lr0,lr)+2*vmd(lr0,lr)  w0=-1  w1=0
    end do
    do lrg=lr+1,norb_dz
      wls = wls+vdint(lr0,lrg)+vdint(lr,lrg)
      wlt = wlt+vdint(lr0,lrg)+vdint(lr,lrg)
      if (jb_sys > 0) then
        wls1 = wls1+vdint(lr0,lrg)+vdint(lr,lrg)
      end if
      !if (jb_sys > 1) then
      !  wlt1 = wlt1+vdint(lr0,lrg)+vdint(lr,lrg)
      !end if
    end do
    do ipae_=1,25
      ipae = ipae_ ! ipae is in common block, is this necessary?
      if (nu_ae(ipae) == 0) cycle
      iwdownv = iw_downwei(jpat,ipae)
      do iwa=0,iwdownv-1
        iwad = iwalk_ad(jpat,ipae,iwa,iwt)
        call prodel(1,wlt,0,ipae,iwad)
      end do
      iwdownv = iw_downwei(jpas,ipae)
      do iwa=0,iwdownv-1
        iwad = iwalk_ad(jpas,ipae,iwa,iws)
        call prodel(1,wls,0,ipae,iwad)
      end do
      if (jb_sys > 0) then
        iwdownv = iw_downwei(jpas,ipae)
        do iwa=0,iwdownv-1
          iwad = iwalk_ad(jpas,ipae,iwa,iws1)     ! s 1-2
          call prodel(1,wls1,0,ipae,iwad)
        end do
      end if
      if (jb_sys > 1) then
        iwdownv = iw_downwei(jpat1,ipae)
        do iwa=0,iwdownv-1
          iwad = iwalk_ad(jpat1,ipae,iwa,iwt)     ! t 1-1
          call prodel(1,wlt,0,ipae,iwad)
        end do
      end if
    end do
  end do
end do

! ------------- end of delm --------------------------------------------
return

end subroutine diagonal_dbl

subroutine diagonal_ext()

#include "drt_h.fh"
#include "intsort_h.fh"
#include "pl_structure_h.fh"
!common /casrst/ja(max_node),jb(max_node),jm(0:max_node),jj(4,0:max_node),kk(0:max_node),no(0:max_innorb),jv,jd(8),jt(8),js(8)
!common /sub_drt/jpad,jpae,ipae,ndim,nohy,ihy(max_wei),jj_sub(4,0:max_node),iy(4,0:max_node),jphy(max_node)

jws0 = 0
do mra=1,8
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
    ! d(1)_800
    lra = norb_all-la+1
    wld = voint(lra,lra)
    call prodel(2,wld,0,ipae,jw)
  end do
end do

!zz = ' out_800_s'
!jps = js(1)

jweis = jws0
do la=1,norb_ext
  !jpd = jd(mra)
  lra = norb_all-la+1
  jweis = jweis+1
  ! (3)_800
  wls = 2.d0*voint(lra,lra)+vdint(lra,lra)
  call prodel(2,wls,0,18,jweis)
end do
do im=1,8
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
      mr = mul_tab(ima,imb)
      if (mr /= im) goto 600
      !jps = js(mr)
      !jpt = jt(mr)
      jws = jws+1
      jwt = jwt+1
      wls = voint(lra,lra)+voint(lrb,lrb)
      wlt = wls
      wls = wls+voint(lrb,lra)+vdint(lrb,lra)   ! w0=-1/2  w1=-3/2
      wlt = wlt-voint(lrb,lra)+vdint(lrb,lra)   ! w0=-1/2  w1=1/2

      call prodel(2,wls,0,ipas,jws)
      call prodel(2,wlt,0,ipat,jwt)
600   continue
      continue
    end do
  end do
end do

! ------------- end of h_delm ------------------------------------------
return

end subroutine diagonal_ext

subroutine prodel(idb,wl,mg1,mg2,mg3)

#include "drt_h.fh"

if (log_prod == 3) then
  call prodel_pt(idb,wl,mg1,mg2,mg3) ! pt1
else
  call prodel_hd(idb,wl,mg1,mg2,mg3)! form_h
end if

return

end subroutine prodel

!  2001.10.2 for norb_act<>0                          mg1,mg2,mg3:
!  idb=1  in dbl_space         ity_up=0-5             jd_type,jd_im,iwd
!  idb=2  in ext_space         ity_down=0-3           je_type,je_im,iwe
!  idb=3  in act_space         ity_up=0-5,itdown=0,3  jp,     mpe,  iwa
!  idb=4  between dbl and act  ity_up=0-5,itdown=0,3  mpe,    iwd,  iwa
!  idb=5  between act and ext  ity_down=0-3           jp,     iwa,  iwe
!  idb=6  between dbl and ext  ity_down=0-3           iwd,    iwa,  iwe
subroutine prodel_hd(idb,wl,mg1,mg2,mg3)

#include "drt_h.fh"
#include "intsort_h.fh"
#include "pl_structure_h.fh"
!common /sub_drt/jpad,jpae,ipae,ndim,nohy,ihy(max_wei),jj_sub(4,0:max_node),iy(4,0:max_node),jphy(max_node)

!ndr = 1
goto(100,200,300,400,500,600),idb
! in dbl_space
100 continue
ipae = mg2
iwad = mg3
isegdownwei = iseg_downwei(ipae)
do mm=iwad+1,iwad+isegdownwei
  vector1(mm) = vector1(mm)+wl
  !if (mm == ndr) then
  !  write(6,'(a8,3i6,2f20.14)') ' in dbl _',mg1,mg2,mg3,wl,vector1(mm)
  !end if
end do
goto 1000
! in ext_space
200 continue
ipae = mg2
iwe = mg3
do jdbl=1,mxnode
  if (nu_ad(jdbl) == 0) cycle
  iw = iw_downwei(jdbl,ipae)
  iwupwei = jpad_upwei(jdbl)
  do iwa=0,iw-1
    do iwd=0,iwupwei-1
      mm = iwalk_ad(jdbl,ipae,iwa,iwd)+iwe
      vector1(mm) = vector1(mm)+wl
      !if(mm == ndr) then
      !  write(nf2,'(a8,3i6,2f20.14)') ' in ext _',mg1,mg2,mg3,wl,vector1(mm)
      !  write(6,'(a8,3i6,2f20.14)') ' in ext _',mg1,mg2,mg3,wl,vector1(mm)
      !end if
    end do
  end do
end do
goto 1000
! in act_space
300 continue
jp = mg1
mpe = mg2
jw = mg3
iwupwei = jpad_upwei(jpad)
isegdownwei = iseg_downwei(ipae)
jph = jphy(jp)
in = ihy(jph)
lwnu = iy(1,mpe)
do jwu=jph+1,jph+in
  iwa = jw+ihy(jwu)-1
  do jwd=1,lwnu
    iwa = iwa+1
    do iwd=0,iwupwei-1
      iwad = iwalk_ad(jpad,ipae,iwa,iwd)
      do iwe=1,isegdownwei
        mm = iwe+iwad
        vector1(mm) = vector1(mm)+wl
        !if (mm == ndr) then
        !  write(6,'(a8,3i6,2f20.14)') ' in act _',mg1,mg2,mg3,wl,vector1(mm)
        !end if
      end do
    end do
  end do
end do
goto 1000
! between dbl and act
400 continue
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
    vector1(mm) = vector1(mm)+wl
    !if (mm == ndr) then
    !  write(6,'(a8,3i6,2f20.14)') ' dbl_act ',mg1,mg2,mg3,wl,vector1(mm)
    !end if
  end do
end do
goto 1000
! between act and ext
500 continue
jp = mg1
iwa0 = mg2
iwe = mg3
iwupwei = jpad_upwei(jpad)
jph = jphy(jp)
in = ihy(jph)
do jwu=jph+1,jph+in
  iwa = iwa0+ihy(jwu)
  do iwd=0,iwupwei-1
    iwad = iwalk_ad(jpad,ipae,iwa,iwd)
    mm = iwe+iwad
    vector1(mm) = vector1(mm)+wl
    !if (mm == ndr) then
    !  write(nf2,'(a8,3i6,2f20.14)') ' act_ext ',mg1,mg2,mg3,wl,vector1(mm)
    !  write(6,'(a8,3i6,2f20.14)') ' act_ext ',mg1,mg2,mg3,wl,vector1(mm)
    !end if
  end do
end do
goto 1000
! between dbl and ext
600 continue
iwd = mg1
iwa = mg2
iwe = mg3
iwad = iwalk_ad(jpad,ipae,iwa,iwd)   ! between dbl,act and ext
mm = iwe+iwad
vector1(mm) = vector1(mm)+wl
!if (mm == ndr) then
!  write(nf2,'(a8,3i6,2f20.14)') ' dbl_ext ',mg1,mg2,mg3,wl,vector1(mm)
!  write(6,'(a8,3i6,2f20.14)') ' dbl_ext ',mg1,mg2,mg3,wl,vector1(mm)
!end if

1000 continue
return

end subroutine prodel_hd

subroutine prodel_pt(idb,wl,mg1,mg2,mg3)

#include "drt_h.fh"
#include "intsort_h.fh"
#include "pl_structure_h.fh"
!common /sub_drt/jpad,jpae,ipae,ndim,nohy,ihy(max_wei),jj_sub(4,0:max_node),iy(4,0:max_node),jphy(max_node)

!if ((jpad == 18) .and. (ipae == 2)) ndr=ndim_h0+1

!ndr = 32257+ndim_h0
goto(100,200,300,400,500,600),idb
! in dbl_space
100 continue
ipae = mg2
iwad = mg3
isegdownwei = jpae_downwei(ipae)
do mm=iwad+1,iwad+isegdownwei
  vector1(mm) = vector1(mm)+wl
  !if (mm == ndr) then
  !  write(nf2,'(a8,3i6,2f20.14)') ' in dbl _',mg1,mg2,mg3,wl,vector1(mm)
  !  write(6,'(a8,3i6,2f20.14)') ' in dbl _',mg1,mg2,mg3,wl,vector1(mm)
  !end if
end do
goto 1000
! in ext_space
200 continue
ipae = mg2
iwe = mg3
iw = ndim
!iw = iseg_dim(jpad,ipae)
iwupwei = jpad_upwei(jpad)
do iwa=0,iw-1
  do iwd=0,iwupwei-1
    mm = iwalk_ad(jpad,ipae,iwa,iwd)+iwe
    vector1(mm) = vector1(mm)+wl
    !if (mm == ndr) then
    !  write(nf2,'(a8,3i6,2f20.14)') ' in ext _',mg1,mg2,mg3,wl,vector1(mm)
    !  write(6,'(a8,3i6,2f20.14)') ' in ext _',mg1,mg2,mg3,wl,vector1(mm)
    !end if
  end do
end do
goto 1000
! in act_space
300 continue
jp = mg1
mpe = mg2
jw = mg3
iwupwei = jpad_upwei(jpad)
isegdownwei = jpae_downwei(ipae)
jph = jphy(jp)
in = ihy(jph)
lwnu = iy(1,mpe)
do jwu=jph+1,jph+in
  iwa = jw+ihy(jwu)-1
  do jwd=1,lwnu
    iwa = iwa+1
    do iwd=0,iwupwei-1
      iwad = iwalk_ad(jpad,ipae,iwa,iwd)
      do iwe=1,isegdownwei
        mm = iwe+iwad
        vector1(mm) = vector1(mm)+wl
        !if (mm == ndr) then
        !  write(nf2,'(a8,3i6,2f20.14)') ' in act _',mg1,mg2,mg3,wl,vector1(mm)
        !  write(6,'(a8,3i6,2f20.14)') ' in act _',mg1,mg2,mg3,wl,vector1(mm)
        !end if
      end do
    end do
  end do
end do
goto 1000
! between dbl and act
400 continue
mpe = mg1
iwd = mg2
iwa = mg3-1
isegdownwei = jpae_downwei(ipae)       ! between dbl and act
jwnu = iy(1,mpe)
do ii=1,jwnu
  iwa = iwa+1
  iwad = iwalk_ad(jpad,ipae,iwa,iwd)
  do iwe=1,isegdownwei
    mm = iwe+iwad                  ! iwl=iwalk_ad
    vector1(mm) = vector1(mm)+wl
    !if (mm == ndr) then
    !  write(nf2,'(a8,3i6,2f20.14)') ' dbl_act ',mg1,mg2,mg3,wl,vector1(mm)
    !  write(6,'(a8,3i6,2f20.14)') ' dbl_act ',mg1,mg2,mg3,wl,vector1(mm)
    !end if
  end do
end do
goto 1000
! between act and ext
500 continue
jp = mg1
iwa0 = mg2
iwe = mg3
iwupwei = jpad_upwei(jpad)
jph = jphy(jp)
in = ihy(jph)
do jwu=jph+1,jph+in
  iwa = iwa0+ihy(jwu)
  do iwd=0,iwupwei-1
    iwad = iwalk_ad(jpad,ipae,iwa,iwd)
    mm = iwe+iwad
    vector1(mm) = vector1(mm)+wl
    !if (mm == ndr) then
    !  write(nf2,'(a8,3i6,2f20.14)') ' act_ext ',mg1,mg2,mg3,wl,vector1(mm)
    !  write(6,'(a8,3i6,2f20.14)') ' act_ext ',mg1,mg2,mg3,wl,vector1(mm)
    !end if
  end do
end do
goto 1000
! between dbl and ext
600 continue
iwd = mg1
iwa = mg2
iwe = mg3
iwad = iwalk_ad(jpad,ipae,iwa,iwd)       ! between dbl,act and ext
mm = iwe+iwad
vector1(mm) = vector1(mm)+wl
!if (mm == ndr) then
!  write(nf2,'(a8,3i6,2f20.14)') ' dbl_ext ',mg1,mg2,mg3,wl,vector1(mm)
!  write(6,'(a8,3i6,2f20.14)') ' dbl_ext ',mg1,mg2,mg3,wl,vector1(mm)
!end if
1000 continue

return

end subroutine prodel_pt

subroutine get_jp(ity,nms,jp,id)

#include "drt_h.fh"

ms = nms
if (id == 1) ms = mul_tab(nms,ns_sm)
goto(10,20,30,40,50,60),ity
10 continue
jp = 1
return
20 continue
jp = 1+ms
return
30 continue
jp = 9+ms
return
40 continue
jp = 17+ms
return
50 continue
jp = 25+ms
return
60 continue
jp = 33+ms

return

end subroutine get_jp
