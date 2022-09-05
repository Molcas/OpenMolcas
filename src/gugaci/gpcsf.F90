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

! generate and print csfs
subroutine found_a_config(ndl,de,npr)

use gugaci_global, only: ipae, iseg_downwei, iw_downwei, iw_sta, jd, jpad, jpad_upwei, jpae, js, jt, jv, map_orb_order, mxnode, &
                         ndim, ndr, ng_sm, nlsm_frz, noidx, norb_act, norb_all, norb_dz, norb_ext, nu_ad, nu_ae, nwalk
                         !, lsmorb, norb_frz
use stdalloc, only: mma_allocate, mma_deallocate
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: ndl, npr
real(kind=wp), intent(in) :: de
integer(kind=iwp) :: i, im, iwupwei, j, jaedownwei, jpad_, l, ndimsum, ne, no_dz, noi, ns, nst
real(kind=wp) :: sqde
integer(kind=iwp), allocatable :: norbindex(:), norbsymm(:), norbwalk(:), nwalk_gamess(:)
character(len=18) :: form1

ndr = ndl
nst = norb_all
do l=1,norb_ext+norb_act
  nwalk(l) = 0
end do
no_dz = nst-norb_dz+1
do l=no_dz,nst
  nwalk(l) = 3
end do
ndimsum = 0
if ((npr == 1) .or. (npr == 0)) then
  jpae = jv
  ipae = 1
  jaedownwei = iseg_downwei(ipae)
  jpad = 1
  iw_sta(jpad,ipae) = ndimsum
  call seg_drt()                             !   jpad_upwei(*)        = jp
  iwupwei = jpad_upwei(jpad)                 !   iw_sta(jpad,jpae)
  iw_downwei(jpad,ipae) = ndim               !   iw_downwei(jpad,jpae)
  ndimsum = ndimsum+ndim*jaedownwei*iwupwei  !   iseg_dim(jpae)    =
  call config_act()                          !   jpae_upwei(jpae)  =
else

  jpae = jv
  ipae = 1
  jaedownwei = iseg_downwei(ipae)
  do jpad_=1,mxnode
    jpad = jpad_ ! jpad is in global module, is this necessary?
    iw_sta(jpad,ipae) = ndimsum
    if (nu_ad(jpad) == 0) cycle
    call seg_drt()                             !   jpad_upwei(*)        = jp
    iwupwei = jpad_upwei(jpad)                 !   iw_sta(jpad,jpae)
    iw_downwei(jpad,ipae) = ndim               !   iw_downwei(jpad,jpae)
    ndimsum = ndimsum+ndim*jaedownwei*iwupwei  !   iseg_dim(jpae)    =
    if (ndim == 0) cycle                       !   jpae_sta(jpae)    = jpa
    call config_act()                          !   jpae_upwei(jpae)  =
  end do
  do im=1,8
    jpae = jd(im)
    ipae = 1+im
    iw_sta(1:mxnode,ipae) = ndimsum
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
      call config_act()
    end do
  end do
  do im=1,8
    jpae = jt(im)
    ipae = 9+im
    iw_sta(1:mxnode,ipae) = ndimsum
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
      call config_act()
    end do
  end do
  do im=1,8
    jpae = js(im)
    ipae = 17+im
    if (nu_ae(ipae) == 0) cycle
    do jpad_=1,mxnode
      jpad = jpad_ ! jpad is in global module, is this necessary?
      if (nu_ad(jpad) == 0) cycle
      call seg_drt()
      if (ndim == 0) cycle
      call config_act()
    end do
  end do
  call config_dbl()
  call config_ext()

end if

if (npr == 0) return

call mma_allocate(norbindex,norb_all,label='norbindex')
call mma_allocate(norbsymm,norb_all,label='norbsymm')
call mma_allocate(norbwalk,norb_all,label='norbwalk')
call mma_allocate(nwalk_gamess,norb_all,label='nwalk_gamess')

nwalk_gamess(1:norb_all) = 0
do i=1,nst
  nwalk_gamess(i) = nwalk(nst-i+1)
end do
if (npr == 1) then
  write(u6,1000) ndr,de
else
  sqde = de*de
  write(u6,1001) ndr,de,sqde
end if

!if (intgen == 1) then
!  ! gamess integral
!  i = 0
!  j = 0
!  do i=norb_frz+1,nst
!    if (nwalk_gamess(map_orb_order(i)) > 0) then
!      j = j+1
!      norbindex(j) = i
!      norbsymm(j) = lsmorb(i)
!      norbwalk(j) = nwalk_gamess(map_orb_order(i))
!    end if
!  end do
!
!  !write(form1,2010) '(4x,i4,1x,',ng_sm,'(1x','i'
!  write(form1,1010) '(4x,a4,',j,'(1x,i3))'
!  write(u6,form1) 'norb',(norbindex(i),i=1,j)
!  write(u6,form1) 'sym ',(norbsymm(i),i=1,j)
!  write(u6,form1) 'walk',(norbwalk(i),i=1,j)
!else
i = 0
j = 0
do im=1,ng_sm
  ns = noidx(im)+nlsm_frz(im)+1
  if (im == ng_sm) then
    ne = norb_all
  else
    ne = noidx(im+1)
  end if
  do i=ns,ne
    if (nwalk_gamess(map_orb_order(i)) > 0) then
      j = j+1
      noi = i-noidx(im)
      norbindex(j) = noi
      norbsymm(j) = im
      norbwalk(j) = nwalk_gamess(map_orb_order(i))
    end if
  end do
end do

write(form1,1010) '(4x,a4,',j,'(1x,i3))'
write(u6,form1) 'norb',(norbindex(i),i=1,j)
write(u6,form1) 'sym ',(norbsymm(i),i=1,j)
write(u6,form1) 'walk',(norbwalk(i),i=1,j)
!end if
call mma_deallocate(norbindex)
call mma_deallocate(norbsymm)
call mma_deallocate(norbwalk)
call mma_deallocate(nwalk_gamess)

return

1000 format(/4x,'csf',i8,6x,'scf energy',f18.8,/)
1001 format(/4x,'csf',i8,2x,'coef',f10.6,2x,'weight',1x,f8.6/)
1010 format(a7,i3,a8)

end subroutine found_a_config

subroutine config_act()

use gugaci_global, only: iy, jj_sub, no, norb_dz, norb_inn
use Definitions, only: iwp

implicit none
integer(kind=iwp) :: idl, jp, jp0, jp1, jw, lr, lr0, mpe

!write(u6,*) '               ***** start h-diaelm *****'
!write(u6,*) '   diagonal_act_d:',jpad,ipae
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
      jw = iy(idl,jp)
      !call prodel(3,wt,jp,mpe,jw)
      call prodel_conf(3,jp,mpe,jw,lr0,0,idl-1)
    end do
    mpe = jj_sub(4,jp)
    if (mpe /= 0) then
      !wt = vdint(lr0,lr0)+Two*voint(lr0,lr0)     !idl=4 hnil=2
      jw = iy(4,jp)
      !call prodel(3,wt,jp,mpe,jw)
      call prodel_conf(3,jp,mpe,jw,lr0,0,3)
    end if
  end do
end do

return

end subroutine config_act

subroutine config_dbl()

use gugaci_global, only: ipae, iw_downwei, jb_sys, jpad, jud, just, lsm_inn, norb_dbl, norb_dz, norb_frz, ns_sm, nu_ae
use Symmetry_Info, only: Mul
use Definitions, only: iwp

implicit none
integer(kind=iwp) :: ipae_, iwa, iwad, iwd, iwdownv, iws, iws1, iwt, jpad1, jpas, jpat, jpat1, lr, lr0, mr, mr0
integer(kind=iwp), external :: iwalk_ad

if (norb_dbl == 0) return
do ipae_=1,25
  ipae = ipae_ ! ipae is in global module, is this necessary?
  if (nu_ae(ipae) == 0) cycle
  iwdownv = iw_downwei(1,ipae)
  do iwa=0,iwdownv-1
    !zz = ' doub_800_v'
    iwad = iwalk_ad(1,ipae,iwa,0)
    !call prodel(1,wt0,0,ipae,iwad)
    call prodel_conf(1,0,ipae,iwad,0,0,1)
  end do
end do
!jps = js(1)
do lr0=norb_frz+1,norb_dz
  mr0 = Mul(lsm_inn(lr0),ns_sm)
  iwd = jud(lr0)
  jpad = 1+mr0
  jpad1 = jpad+24
  !wld = wt0-voint(lr0,lr0)-vdint(lr0,lr0)
  !wls = wld-voint(lr0,lr0)
  do ipae_=1,25
    ipae = ipae_ ! ipae is in global module, is this necessary?
    if (nu_ae(ipae) == 0) cycle
    iwdownv = iw_downwei(jpad,ipae)
    do iwa=0,iwdownv-1
      iwad = iwalk_ad(jpad,ipae,iwa,iwd)
      !call prodel(1,wld,0,ipae,iwad)
      call prodel_conf(1,0,ipae,iwad,lr0,0,2)     !d
    end do
  end do

  if (jb_sys > 0) then
    do ipae_=1,25
      ipae = ipae_ ! ipae is in global module, is this necessary?
      if (nu_ae(ipae) == 0) cycle
      iwdownv = iw_downwei(jpad1,ipae)
      do iwa=0,iwdownv-1
        iwad = iwalk_ad(jpad1,ipae,iwa,iwd)
        !call prodel(1,wld,0,ipae,iwad)
        call prodel_conf(1,0,ipae,iwad,lr0,0,5)     !dd
      end do
    end do
  end if
  jpad = 17+ns_sm
  iwd = just(lr0,lr0)
  do ipae_=1,25
    ipae = ipae_ ! ipae is in global module, is this necessary?
    if (nu_ae(ipae) == 0) cycle
    iwdownv = iw_downwei(jpad,ipae)
    do iwa=0,iwdownv-1
      iwad = iwalk_ad(jpad,ipae,iwa,iwd)
      !call prodel(1,wls,0,ipae,iwad)
      call prodel_conf(1,0,ipae,iwad,lr0,0,4)     !s
    end do
  end do

  if (lr0 == norb_dz) cycle
  !wld0 = wld

  do lr=lr0+1,norb_dz
    mr = Mul(mr0,lsm_inn(lr))
    jpat = 9+mr
    jpas = 17+mr
    jpat1 = jpat+24
    iws = just(lr0,lr)
    iwt = iws
    !wld = wld0-voint(lr,lr)-vdint(lr,lr)
    do ipae_=1,25
      ipae = ipae_ ! ipae is in global module, is this necessary?
      if (nu_ae(ipae) == 0) cycle
      iwdownv = iw_downwei(jpat,ipae)
      do iwa=0,iwdownv-1
        iwad = iwalk_ad(jpat,ipae,iwa,iwt)
        !call prodel(1,wld,0,ipae,iwad)
        call prodel_conf(1,0,ipae,iwad,lr0,lr,3)   !t
      end do
      if (jb_sys > 1) then
        iwdownv = iw_downwei(jpat1,ipae)
        do iwa=0,iwdownv-1
          iwad = iwalk_ad(jpat1,ipae,iwa,iwt)         !tt
          !call prodel(1,wld,0,ipae,iwad)
          call prodel_conf(1,0,ipae,iwad,lr0,lr,6)
        end do
      end if
      iwdownv = iw_downwei(jpas,ipae)
      do iwa=0,iwdownv-1
        iwad = iwalk_ad(jpas,ipae,iwa,iws)
        !call prodel(1,wld,0,ipae,iwad)
        call prodel_conf(1,0,ipae,iwad,lr0,lr,4)   !s
      end do
      if (jb_sys > 0) then
        iws1 = just(lr,lr0)
        do iwa=0,iwdownv-1
          iwad = iwalk_ad(jpas,ipae,iwa,iws1)
          !call prodel(1,wld,0,ipae,iwad)
          call prodel_conf(1,0,ipae,iwad,lr0,lr,7)      !ss
        end do
      end if
    end do
  end do
end do

return

end subroutine config_dbl

subroutine config_ext()

use gugaci_global, only: ibsm_ext, iesm_ext, ipae, lsm, norb_ext
use Symmetry_Info, only: Mul
use Definitions, only: iwp

implicit none
integer(kind=iwp) :: im, ima, imb, ipas, ipat, jw, jweis, jws, jws0, jwt, la, laend, lasta, lb, lrzz, mr, mra

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
    !lra = norb_all-la+1
    !wld = voint(lra,lra)
    !call prodel(2,wld,0,ipae,jw)
    call prodel_conf(2,0,ipae,jw,la,0,2)      !d
  end do
end do

!zz = ' out_800_s'
!jps = js(1)

jweis = jws0
do la=1,norb_ext
  !jpd = jd(mra)
  !lra = norb_all-la+1
  jweis = jweis+1
  !wls = Two*voint(lra,lra)+vdint(lra,lra)
  !call prodel(2,wls,0,18,jweis)
  call prodel_conf(2,0,18,jweis,la,0,4)      !s
end do
do im=1,8
  jws = 0
  jwt = 0
  ipat = 9+im
  ipas = 17+im
  do la=2,norb_ext
    !lra = norb_all-la+1
    ima = lsm(la)
    do lb=1,la-1
      !lrbi = norb_all-lb+1
      imb = lsm(lb)
      mr = Mul(ima,imb)
      if (mr /= im) cycle
      !jps = js(mr)
      !jpt = jt(mr)
      jws = jws+1
      jwt = jwt+1
      !wls = voint(lra,lra)+voint(lrb,lrb)
      !wlt = wls
      !wls = wls+voint(lrb,lra)+vdint(lrb,lra)   ! w0=-1/2  w1=-3/2
      !wlt = wlt-voint(lrb,lra)+vdint(lrb,lra)   ! w0=-1/2  w1=1/2

      !call prodel(2,wls,0,ipas,jws)
      !call prodel(2,wlt,0,ipat,jwt)
      call prodel_conf(2,0,ipas,jws,la,lb,4)      !s
      call prodel_conf(2,0,ipat,jwt,la,lb,3)      !t
    end do
  end do
end do

return

end subroutine config_ext

! idb=1  in dbl_space      ity_up=0-5             jd_type,jd_im,iwd
! idb=2  in ext_space      ity_down=0-3           je_type,je_im,iwe
! idb=3  in act_space      ity_up=0-5,itdown=0,3  jp,     mpe,  iwa
subroutine prodel_conf(idb,mg1,mg2,mg3,lr01,lr02,jpty)

use gugaci_global, only: ihy, ipae, iseg_downwei, iw_downwei, iy, jpad, jpad_upwei, jphy, mxnode, ndr, norb_all, nu_ad, nwalk
use Definitions, only: iwp

implicit none
integer(kind=iwp), intent(in) :: idb, mg1, mg2, mg3, lr01, lr02, jpty
integer(kind=iwp) :: in_, isegdownwei, iw, iwa, iwad, iwd, iwe, iwupwei, jdbl, jp, jph, jw, jwd, jwu, lr1, lr2, lwnu, mm, mpe
integer(kind=iwp), external :: iwalk_ad

select case (idb)
  case default ! (1)
    ! in dbl_space
    ! jpty=  1,  2,  3,  4,  5,  6,  7
    !        v   d   t   s  dd  tt  ss
    ipae = mg2
    iwad = mg3
    lr1 = norb_all-lr01+1
    lr2 = norb_all-lr02+1
    isegdownwei = iseg_downwei(ipae)
    do mm=iwad+1,iwad+isegdownwei
      if (mm == ndr) then
        select case (jpty)
          case default ! (1)
          case (2)
            nwalk(lr1) = 2
          case (3)
            nwalk(lr1) = 2
            nwalk(lr2) = 2
          case (4)
            nwalk(lr1) = 2
            nwalk(lr2) = 1
            if (lr02 == 0) nwalk(lr1) = 0
            if (lr01 == 0) nwalk(lr2) = 0
          case (5)
            nwalk(lr1) = 1
          case (6)
            nwalk(lr1) = 1
            nwalk(lr2) = 1
          case (7)
            nwalk(lr1) = 1
            nwalk(lr2) = 2
        end select
        exit
      end if
    end do

  case (2)
    ! in ext_space
    ! jpty=  1,  2,  3,  4
    !        v   d   t   s
    ipae = mg2
    iwe = mg3
    lr1 = lr01
    lr2 = lr02
    outer1: do jdbl=1,mxnode
      if (nu_ad(jdbl) == 0) cycle
      iw = iw_downwei(jdbl,ipae)
      iwupwei = jpad_upwei(jdbl)
      do iwa=0,iw-1
        do iwd=0,iwupwei-1
          mm = iwalk_ad(jdbl,ipae,iwa,iwd)+iwe
          if (mm == ndr) then
            select case (jpty)
              case default ! (1)
              case (2)
                nwalk(lr1) = 1
              case (3)
                nwalk(lr1) = 1
                nwalk(lr2) = 1
              case (4)
                nwalk(lr1) = 2
                nwalk(lr2) = 1
                if (lr02 == 0) nwalk(lr1) = 3
            end select
            exit outer1
          end if
        end do
      end do
    end do outer1

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
    lr1 = norb_all-lr01+1
    outer2: do jwu=jph+1,jph+in_
      iwa = jw+ihy(jwu)-1
      do jwd=1,lwnu
        iwa = iwa+1
        do iwd=0,iwupwei-1
          iwad = iwalk_ad(jpad,ipae,iwa,iwd)
          do iwe=1,isegdownwei
            mm = iwe+iwad
            if (mm == ndr) then
              nwalk(lr1) = jpty
              exit outer2
            end if
          end do
        end do
      end do
    end do outer2
end select

return

end subroutine prodel_conf
