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

! drt related subroutines

subroutine active_drt()

use gugaci_global, only: iseg_downwei, iseg_sta, iseg_upwei, jd, jj, jpad_upwei, js, jt, jv, LuDrt, max_node, mxnode, nci_dim, &
                         ng_sm, no, norb_act, norb_inn, nu_ad, nu_ae !, logic_mr, logic_mrelcas
use stdalloc, only: mma_allocate, mma_deallocate
use Definitions, only: iwp, u6

implicit none
integer(kind=iwp) :: i, im, iseg_dim(25), jde, jdim, jdn, jds, jji, jp, jpe, jpn, jsim, jtim, ndi
integer(kind=iwp), allocatable :: iin(:)

nci_dim = 0
iseg_dim(:) = 0
if (norb_act == 0) then
  !====================  norb_act=0 ====================================
  iseg_sta(1) = 0
  iseg_dim(1) = 1
  do im=1,ng_sm
    jdim = nu_ae(1+im)
    jtim = nu_ae(9+im)
    jsim = nu_ae(17+im)
    jd(im) = jdim
    jt(im) = jtim
    js(im) = jsim
    iseg_dim(jdim) = iseg_downwei(jdim)*jpad_upwei(jdim)
    iseg_dim(jtim) = iseg_downwei(jtim)*jpad_upwei(jtim)
    iseg_dim(jsim) = iseg_downwei(jsim)*jpad_upwei(jsim)
    if (iseg_dim(jdim) == 0) then
      jd(im) = 0
      nu_ad(jdim) = 0
      nu_ae(jdim) = 0
    end if
    if (iseg_dim(jtim) == 0) then
      jt(im) = 0
      nu_ad(jtim) = 0
      nu_ae(jtim) = 0
    end if
    if (iseg_dim(jsim) == 0) then
      js(im) = 0
      nu_ad(jsim) = 0
      nu_ae(jsim) = 0
    end if
  end do
  do jp=2,25
    iseg_sta(jp) = iseg_sta(jp-1)+iseg_dim(jp-1)
  end do
  nci_dim = iseg_sta(25)+iseg_dim(25)
  iseg_sta(26) = nci_dim
  do jp=1,25
    if (iseg_downwei(jp) == 0) cycle
    iseg_upwei(jp) = iseg_dim(jp)/iseg_downwei(jp)
  end do
else
  !====================  norb_act<>0 ===================================
  !if (logic_mr) call rst()         !npp=2
  !if (logic_mrelcas) call rcas()   !npp=3

  write(u6,*) ' '
  write(u6,*) ' now reading distinct row tableau'
  call readdrt(ludrt)

  nu_ae(1) = jv
  do im=1,ng_sm
    nu_ae(1+im) = jd(im)
    nu_ae(9+im) = jt(im)
    nu_ae(17+im) = js(im)
  end do
  jds = 1
  jde = mxnode
  jpe = no(norb_inn+1)

  jp = jv
  call mma_allocate(iin,[0,max_node],label='iin')
  !iin(1:jpe) = 0
  !iin(0) = 0
  iin(:) = 0
  iin(jp) = 1
  iseg_sta(1) = 0       ! for node_ext
  iseg_dim(1) = 0
  do jpn=jpe,1,-1
    do i=1,4
      jji = jj(i,jpn)
      if (iin(jji) == 0) cycle
      iin(jpn) = iin(jpn)+iin(jji)
    end do
  end do
  do jdn=jds,jde
    if (nu_ad(jdn) == 0) cycle
    ndi = iin(jdn)*iseg_downwei(1)*jpad_upwei(jdn)
    iseg_dim(1) = iseg_dim(1)+ndi
  end do
  do im=1,ng_sm
    jp = jd(im)
    iseg_sta(1+im) = nci_dim       ! for node_d_ext
    iseg_dim(1+im) = 0
    if (jp == 0) cycle
    iin(1:jpe) = 0
    iin(0) = 0
    iin(jp) = 1
    do jpn=jpe,1,-1
      do i=1,4
        jji = jj(i,jpn)
        if (iin(jji) == 0) cycle
        iin(jpn) = iin(jpn)+iin(jji)
      end do
    end do
    do jdn=jds,jde
      if (nu_ad(jdn) == 0) cycle
      ndi = iin(jdn)*iseg_downwei(1+im)*jpad_upwei(jdn)
      iseg_dim(1+im) = iseg_dim(1+im)+ndi
    end do
  end do

  do im=1,ng_sm
    jp = jt(im)
    iseg_sta(9+im) = nci_dim       ! for node_t_ext
    iseg_dim(9+im) = 0
    if (jp == 0) cycle
    iin(1:jpe) = 0
    iin(0) = 0
    iin(jp) = 1
    do jpn=jpe,1,-1
      do i=1,4
        jji = jj(i,jpn)
        if (iin(jji) == 0) cycle
        iin(jpn) = iin(jpn)+iin(jji)
      end do
    end do
    do jdn=jds,jde
      if (nu_ad(jdn) == 0) cycle
      ndi = iin(jdn)*iseg_downwei(9+im)*jpad_upwei(jdn)
      iseg_dim(9+im) = iseg_dim(9+im)+ndi
    end do
  end do
  do im=1,ng_sm
    jp = js(im)
    iseg_sta(17+im) = nci_dim       ! for node_s_ext
    iseg_dim(17+im) = 0
    if (jp == 0) cycle
    iin(1:jpe) = 0
    iin(0) = 0
    iin(jp) = 1
    do jpn=jpe,1,-1
      do i=1,4
        jji = jj(i,jpn)
        if (iin(jji) == 0) cycle
        iin(jpn) = iin(jpn)+iin(jji)
      end do
    end do
    do jdn=jds,jde
      if (nu_ad(jdn) == 0) cycle
      ndi = iin(jdn)*iseg_downwei(17+im)*jpad_upwei(jdn)
      iseg_dim(17+im) = iseg_dim(17+im)+ndi
    end do
  end do
  call mma_deallocate(iin)
  do jp=2,25
    iseg_sta(jp) = iseg_sta(jp-1)+iseg_dim(jp-1)
    !write(nf2,'(3i10)') jp-1,iseg_sta(jp-1),iseg_dim(jp-1)      !to
  end do
  !write(nf2,'(3i10)') 25,iseg_sta(25),iseg_dim(25)          !to de
  nci_dim = iseg_sta(25)+iseg_dim(25)
  iseg_sta(26) = nci_dim
  do jp=1,25
    if (iseg_downwei(jp) == 0) cycle
    iseg_upwei(jp) = iseg_dim(jp)/iseg_downwei(jp)
  end do
  !call get_ivy()
end if

! to the end of dbl,act,ext parts
call dbl_downwalk()
write(u6,*) '  end of drt,nci_dim= ',nci_dim
!write(u6,*) '-----------------------------------------------'
!write(u6,*) '         ci drt-information'
!write(u6,*) '       num.of orbitals:       ',nst
!write(u6,*) '       num.of froz-orbitals:  ',nzst
!write(u6,*) '       num.of active-orbitals:',ndh
!write(u6,*) '       num.of electrons:      ',nel
!write(u6,*) '       multiplicity:          ',mult
!write(u6,*) '       num.of configurations: ',ndime
!write(u6,*) '       symmetry:              ',zm
!write(u6,*) '-----------------------------------------------'

return

end subroutine active_drt

!subroutine rst()
!
!use gugaci_global, only: LuDrt
!use Definitions, only: u6
!
!implicit none
!
!write(u6,*) ' '
!write(u6,*) ' now reading distinct row tableau'
!call readdrt(ludrt)
!!open(21,file='fort.drt',form='unformatted')
!!read(21) id
!!write(u6,*) ' id=',id
!!read(21) ja(1:id),jb(1:id),jm(1:id)
!!read(21) jj(1:4,0:id)
!!read(21) kk(0:id)
!!read(21) no(0:norb_inn+1)
!!read(21) jv,jd(1:8),jt(1:8),js(1:8)
!!close(21)
!
!return
!
!end subroutine rst

!subroutine ref_gfs(nel,ndj,locu,nm)
!
!use gugaci_global, only: lsm_inn, max_ref, norb_dz, norb_inn, nstart_act, spin
!use Symmetry_Info, only: Mul
!use stdalloc, only: mma_allocate, mma_deallocate
!use Definitions, only: iwp, u6
!
!implicit none
!integer(kind=iwp), intent(in) :: nel, nm
!integer(kind=iwp), intent(out) :: ndj, locu(8,max_ref)
!integer(kind=iwp) :: i, l1, l2, l3, l4, l5, l6, l7, l8, ldj, lh, lhe, lhs, lhsm(8), lm, lpsum, m, m1, m2, m3, m4, m5, m6, m7, m8, &
!                     mdj, mys, ne_act, ne_s, nes, npair, nre
!integer(kind=iwp), allocatable :: lscu(:,:)
!
!ne_act = nel-2*norb_dz
!ne_s = nint(spin*2)
!lhs = nstart_act
!lhe = norb_inn
!lhsm(1:8) = 0
!do lh=lhs,lhe
!  lm = lsm_inn(lh)
!  lhsm(lm) = lhsm(lm)+1
!end do
!mdj = 0
!call mma_allocate(lscu,[0,8],[1,max_ref],label='lscu')
!do nes=ne_s,ne_act,2
!  do l1=0,lhsm(1)
!    do l2=0,lhsm(2)
!      do l3=0,lhsm(3)
!        do l4=0,lhsm(4)
!          do l5=0,lhsm(5)
!            do l6=0,lhsm(6)
!              do l7=0,lhsm(7)
!                do l8=0,lhsm(8)
!                  lpsum = l1+l2+l3+l4+l5+l6+l7+l8
!                  if (lpsum /= nes) cycle
!                  mys = 1
!                  if (mod(l1,2) == 1) mys = Mul(mys,1)
!                  if (mod(l2,2) == 1) mys = Mul(mys,2)
!                  if (mod(l3,2) == 1) mys = Mul(mys,3)
!                  if (mod(l4,2) == 1) mys = Mul(mys,4)
!                  if (mod(l5,2) == 1) mys = Mul(mys,5)
!                  if (mod(l6,2) == 1) mys = Mul(mys,6)
!                  if (mod(l7,2) == 1) mys = Mul(mys,7)
!                  if (mod(l8,2) == 1) mys = Mul(mys,8)
!                  if (mys /= nm) cycle
!                  mdj = mdj+1
!                  lscu(0,mdj) = lpsum
!                  lscu(1,mdj) = l1
!                  lscu(2,mdj) = l2
!                  lscu(3,mdj) = l3
!                  lscu(4,mdj) = l4
!                  lscu(5,mdj) = l5
!                  lscu(6,mdj) = l6
!                  lscu(7,mdj) = l7
!                  lscu(8,mdj) = l8
!                end do
!              end do
!            end do
!          end do
!        end do
!      end do
!    end do
!  end do
!end do
!ndj = 0
!do m=1,mdj
!  npair = (ne_act-lscu(0,m))/2
!  do l1=0,lhsm(1)-lscu(1,m)
!    do l2=0,lhsm(2)-lscu(2,m)
!      do l3=0,lhsm(3)-lscu(3,m)
!        do l4=0,lhsm(4)-lscu(4,m)
!          do l5=0,lhsm(5)-lscu(5,m)
!            do l6=0,lhsm(6)-lscu(6,m)
!              do l7=0,lhsm(7)-lscu(7,m)
!                outer: do l8=0,lhsm(8)-lscu(8,m)
!                  lpsum = l1+l2+l3+l4+l5+l6+l7+l8
!                  if (lpsum == npair) then
!                    m1 = l1*2+lscu(1,m)
!                    m2 = l2*2+lscu(2,m)
!                    m3 = l3*2+lscu(3,m)
!                    m4 = l4*2+lscu(4,m)
!                    m5 = l5*2+lscu(5,m)
!                    m6 = l6*2+lscu(6,m)
!                    m7 = l7*2+lscu(7,m)
!                    m8 = l8*2+lscu(8,m)
!                    do ldj=1,ndj
!                      if ((m1 == locu(1,ldj)) .and. (m2 == locu(2,ldj)) .and. (m3 == locu(3,ldj)) .and. (m4 == locu(4,ldj)) .and. &
!                          (m5 == locu(5,ldj)) .and. (m6 == locu(6,ldj)) .and. (m7 == locu(7,ldj)) .and. (m8 == locu(8,ldj))) &
!                        cycle outer
!                    end do
!                    ndj = ndj+1
!                    locu(1,ndj) = m1
!                    locu(2,ndj) = m2
!                    locu(3,ndj) = m3
!                    locu(4,ndj) = m4
!                    locu(5,ndj) = m5
!                    locu(6,ndj) = m6
!                    locu(7,ndj) = m7
!                    locu(8,ndj) = m8
!                  end if
!                end do outer
!              end do
!            end do
!          end do
!        end do
!      end do
!    end do
!  end do
!end do
!call mma_deallocate(lscu)
!
!do nre=1,ndj
!  write(u6,'(5x,i6,8i3)') nre,(locu(i,nre),i=1,8)
!end do
!
!return
!
!end subroutine ref_gfs

!subroutine rcas()
!
!use gugaci_global, only: LuDrt
!use Definitions, only: u6
!
!implicit none
!
!write(u6,*) ' '
!write(u6,*) ' now reading distinct row tableau'
!call readdrt(ludrt)
!
!!write(u6,*) 'bbs debug rcas,kk(27)',kk(27)
!
!!open(21,file='fort.drt',form='unformatted')
!!read(21) id
!!read(21) ja(1:id),jb(1:id),jm(1:id)
!!read(21) jj(1:4,0:id)
!!read(21) kk(0:id)
!!read(21) no(0:norb_inn+1)
!!read(21) jv,jd(1:8),jt(1:8),js(1:8)
!!close(21)
!
!return
!
!end subroutine rcas

!subroutine check_rcas3(jk,ind,inb,ndj,locu)
!
!use gugaci_global, only: ja, jb, max_node
!use stdalloc, only: mma_allocate, mma_deallocate
!use Definitions, only: iwp
!
!implicit none
!integer(kind=iwp), intent(in) :: jk, ind(8,max_node), ndj, locu(8,ndj)
!integer(kind=iwp), intent(out) :: inb
!integer(kind=iwp) :: i, iex, lsym(8), m, nsumel
!integer(kind=iwp), allocatable :: iexcit(:)
!
!inb = 0
!nsumel = 0
!do i=1,8
!  lsym(i) = ind(i,jk)
!  nsumel = nsumel+lsym(i)
!end do
!call mma_allocate(iexcit,ndj,label='iexcit')
!do i=1,ndj
!  iexcit(i) = 0
!  do m=1,8
!    iex = lsym(m)-locu(m,i)
!    if (iex > 0) then
!      iexcit(i) = iexcit(i)+iex
!    end if
!  end do
!end do
!inb = iexcit(1)
!do i=2,ndj
!  inb = min(inb,iexcit(i))
!end do
!call mma_deallocate(iexcit)
!inb = inb+ja(jk)*2+jb(jk)
!
!return
!
!end subroutine check_rcas3

subroutine irfrst()
! ifrno(j)=i
! irfno(i)=j no. i ref is no. j cfs in h0

use gugaci_global, only: ifrno, iref_occ, irf, irfno, n_ref, nci_h0, norb_act, norb_all, norb_dz, nwalk
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: One
use Definitions, only: iwp, u6

implicit none
integer(kind=iwp) :: i, icount, icsfocc, icsfwlk, ii, ij, im, j, ndimh0
integer(kind=iwp), allocatable :: iwalktmp(:)

icsfwlk = 0
ndimh0 = nci_h0 !iw_sta(2,1)
icount = 0
call mma_allocate(iwalktmp,norb_all,label='iwalktmp')
do i=1,n_ref
  outer: do j=1,ndimh0
    call found_a_config(j,One,0)
    do im=1,norb_all
      iwalktmp(im) = nwalk(norb_all-im+1)
    end do

    ij = norb_dz
    do ii=1,norb_act
      ij = ij+1
      icsfwlk = iwalktmp(ij)
      if (icsfwlk == 3) then
        icsfocc = 2
      else if (icsfwlk == 2) then
        icsfocc = 1
      else if (icsfwlk == 1) then
        icsfocc = 1
      else if (icsfwlk == 0) then
        icsfocc = 0
      else
        icsfocc = -1
      end if
      if (icsfocc /= iref_occ(ij,i)) cycle outer
    end do
    icount = icount+1
    irfno(icount) = j
    ifrno(j) = icount
    !write(u6,2000) iwalktmp(norb_dz+1:norb_dz+norb_act)
    !write(u6,*) 'icount',icount,j
  end do outer
end do
call mma_deallocate(iwalktmp)

irf = icount
write(u6,3000) irf

return

!1000 format(1x,'warnning!the selected csf is not in references states')
!2000 format(1x,'the selected csf is :',2x,32(i1))
3000 format(1x,'number of gelfand states in referance space:',1x,i4)
!...end of irfrst

end subroutine irfrst

!subroutine irfrst_bak(iselcsf_occ)
!
!use gugaci_global, only: ifrno, iref_occ, irf, irfno, max_innorb, max_ref, mjn, mroot, n_ref, nci_h0, norb_act, norb_all, norb_dz, &
!                         nwalk
!use stdalloc, only: mma_allocate, mma_deallocate
!use Constants, only: One
!use Definitions, only: iwp, u6
!
!implicit none
!integer(kind=iwp), intent(in) :: iselcsf_occ(max_innorb,max_ref)
!integer(kind=iwp) :: i, icount, icsfocc, icsfwlk, ii, ij, im, io, iocsf, ire, j, ndimh0, nocc
!logical(kind=iwp) :: log_exist
!integer(kind=iwp), allocatable :: iwalktmp(:)
!
!nocc = 0
!do i=1,mroot
!  log_exist = .false.
!  outer1: do ire=1,n_ref
!    ij = norb_dz
!    do io=1,norb_act
!      ij = ij+1
!      if (iselcsf_occ(io,i) == 3) nocc = 2
!      if (iselcsf_occ(io,i) == 2) nocc = 1
!      if (iselcsf_occ(io,i) == 1) nocc = 1
!      if (iselcsf_occ(io,i) == 0) nocc = 0
!      if (nocc /= iref_occ(io+norb_dz,ire)) cycle outer1
!    end do
!    log_exist = .true.
!    if (log_exist) exit
!  end do outer1
!  if (.not. log_exist) then
!    write(u6,1000)
!    write(u6,2000) iselcsf_occ(1:norb_act,i)
!    write(u6,*) ' please select this state as reference state'
!  end if
!end do
!icsfocc = 0
!ndimh0 = nci_h0 !iw_sta(2,1)
!icount = 0
!call mma_allocate(iwalktmp,norb_all,label='iwalktmp')
!do i=1,n_ref
!  outer2: do j=1,ndimh0
!    call found_a_config(j,One,0)
!    do im=1,norb_all
!      iwalktmp(im) = nwalk(norb_all-im+1)
!    end do
!
!    log_exist = .false.
!    ij = norb_dz
!    do ii=1,norb_act
!      ij = ij+1
!      icsfwlk = iwalktmp(ij)
!      if (icsfwlk == 3) icsfocc = 2
!      if (icsfwlk == 2) icsfocc = 1
!      if (icsfwlk == 1) icsfocc = 1
!      if (icsfwlk == 0) icsfocc = 0
!      if (icsfocc /= iref_occ(ij,i)) cycle outer2
!    end do
!    log_exist = .true.
!    icount = icount+1
!    irfno(icount) = j
!    ifrno(j) = icount
!    !write(u6,2000) iwalktmp(norb_dz+1:norb_dz+norb_act)
!    !write(u6,*) 'icount',icount,j
!  end do outer2
!end do
!call mma_deallocate(iwalktmp)
!
!irf = icount
!do i=1,2*mroot
!  iocsf = mjn(i)
!  log_exist = .false.
!  do j=1,icount
!    if (iocsf == irfno(j)) then
!      log_exist = .true.
!      exit
!    end if
!  end do
!  if (.not. log_exist) then
!    irf = irf+1
!    irfno(irf) = iocsf
!    ifrno(iocsf) = irf
!  end if
!end do
!
!write(u6,3000) irf
!
!return
!
!1000 format(1x,'warnning!the selected csf is not in references states')
!2000 format(1x,'the selected csf is :',2x,32(i1))
!3000 format(1x,'number of gelfand states in referance space:',1x,i4)
!!...end of irfrst_bak
!
!end subroutine irfrst_bak

!function min_itexcit(indjk)
!
!use gugaci_global, only: ndjgrop, ndjmod
!use Definitions, only: iwp
!
!implicit none
!integer(kind=iwp) :: min_itexcit
!integer(kind=iwp), intent(in) :: indjk(4)
!integer(kind=iwp) :: indexcit, ixcit, lref, ngrop, nj
!! integer*4 indjk  =  00 00 00 00 00 00 00 00 00 00  00 00 00 00 00
!! indexcit=  ir1 ir2 ir3 ir4 ir5 ir6 ir7 ir8 ......... ir15
!
!!-----------------------------------------------------------------------
!min_itexcit = 3
!nj = 0
!do ngrop=1,ndjgrop-1         !1-(ndjgrop-1) grop
!  indexcit = indjk(ngrop)
!  do lref=1,15
!    ixcit = ishft(indexcit,-2*(lref-1))
!    ixcit = mod(ixcit,4)
!    !if (ixcit /= 0) ixcit = ixcit-1
!    !if (ixcit == 0) ixcit = 3
!    min_itexcit = min(min_itexcit,ixcit)
!    if (min_itexcit == 0) return
!  end do
!  nj = nj+15
!end do
!indexcit = indjk(ndjgrop)       !last grop
!do lref=1,ndjmod
!  ixcit = ishft(indexcit,-2*(lref-1))
!  ixcit = mod(ixcit,4)
!  !if (ixcit /= 0) ixcit = ixcit-1
!  !if (ixcit == 0) ixcit = 3
!  min_itexcit = min(min_itexcit,ixcit)
!  if (min_itexcit == 0) return
!end do
!!-----------------------------------------------------------------------
!return
!
!end function min_itexcit

!subroutine njexcit(idcc,indjk,locuk0,n_ref)
!
!use gugaci_global, only: ndjgrop, ndjmod
!use Definitions, only: iwp
!
!implicit none
!integer(kind=iwp), intent(in) :: idcc, n_ref, locuk0(n_ref)
!integer(kind=iwp), intent(inout) :: indjk(4)
!integer(kind=iwp) :: indexcit, ixcit, lref, ngrop, nj
!
!nj = 0
!do ngrop=1,ndjgrop-1         !1-(ndjgrop-1) grop
!  indexcit = indjk(ngrop)
!  indjk(ngrop) = 0
!  do lref=1,15
!    ixcit = ishft(indexcit,-2*(lref-1))
!    ixcit = mod(ixcit,4)
!    if (idcc-locuk0(nj+lref) == 1) ixcit = ixcit+1
!    if (idcc-locuk0(nj+lref) == 2) ixcit = ixcit+2
!    if (ixcit >= 3) ixcit = 3
!    indjk(ngrop) = ishft(ixcit,2*(lref-1))+indjk(ngrop)
!  end do
!  nj = nj+15
!end do
!indexcit = indjk(ndjgrop)       !last grop
!indjk(ndjgrop) = 0
!do lref=1,ndjmod
!  ixcit = ishft(indexcit,-2*(lref-1))
!  ixcit = mod(ixcit,4)
!  if (idcc-locuk0(nj+lref) == 1) ixcit = ixcit+1
!  if (idcc-locuk0(nj+lref) == 2) ixcit = ixcit+2
!  if (ixcit >= 3) ixcit = 3
!  indjk(ngrop) = ishft(ixcit,2*(lref-1))+indjk(ngrop)
!end do
!
!return
!
!end subroutine njexcit
