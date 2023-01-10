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

subroutine geth0()

use gugaci_global, only: ecih0, escf, indx, ipae, irf, irfno, jpad, jpadl, jpae, logic_mr, max_h0, max_innorb, max_kspace, &
                         max_ref, max_root, mroot, ndim, ndim_h0, norb_act, nu_ad, nu_ae, vcm, vd, ve, vector1, vector2, vu
                         !, len_str, tmpdir
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: One
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp) :: i, ijm, kval, l, m, mmspace, mn, ndim0, nxb, nxh
real(kind=wp) :: vad0(max_h0)
!character(len=256) :: filename
integer(kind=iwp), allocatable :: iselcsf_occ(:,:)
real(kind=wp), allocatable :: vb1(:), vb2(:)
real(kind=wp), parameter :: dcrita = 5.0e-6_wp

call mma_allocate(iselcsf_occ,max_innorb,max_ref,label='iselcsf_occ')
if (.not. logic_mr) then
  call minevalue(iselcsf_occ)
end if
!=======================================================================
! calculate ndim_h0
if (norb_act == 0) then
  ndim_h0 = 1
  vector2(1) = vector1(1)
  irfno(1) = 1
  !return
else
  ipae = 1
  jpae = nu_ae(ipae)
  if (jpae == 0) return
  jpadl = 1
  if (nu_ad(jpadl) == 0) return
  jpad = jpadl
  call seg_drt()
  if (ndim == 0) return
  ndim_h0 = ndim

  if (mroot > ndim_h0) then
    write(u6,*) '    mroot> ndim_h0, mroot,ndim_h0=',mroot,ndim_h0
    mroot = min(ndim_h0,mroot)
    write(u6,*) '   ',mroot,'roots are calculated'
  end if

  call copy_to_drtl()

  if (logic_mr) then
    call irfrst()
    if (mroot > ndim_h0) then
      write(u6,*) '    mroot> ndim_h0, mroot,ndim_h0=',mroot,ndim_h0
      mroot = min(ndim_h0,mroot)
      write(u6,*) '   ',mroot,'roots are calculated'
    end if
    call minevalue(iselcsf_occ)
  end if
end if
call mma_deallocate(iselcsf_occ)

!=======================================================================
!if (logic_mr) ndim_h0 = irf
if (.not. logic_mr) then
  ndim0 = ndim_h0
else
  ndim0 = irf
end if
write(u6,*) '     ================================'
write(u6,*) '         step 1: diagnalization h0   '
write(u6,*) '            ndim_h0=',ndim0
write(u6,*) '     ================================'

if (ndim_h0 == 1) then
  ecih0(1) = escf(1)
  vcm(1) = One
else

  call formh0()   ! for log_mr, ndim_h0 changed to irf in this subroutine
  !=====================================================================
  if (allocated(vcm)) call mma_deallocate(vcm)
  call mma_allocate(vcm,ndim_h0*mroot,label='vcm')

  if (ndim_h0 <= 30) then
    call hotred(max_kspace,ndim_h0,vector2,vd,ve,vu)
    call qlcm(max_kspace,ndim_h0,vd,ve,vu)
    ijm = 0
    do m=1,mroot
      ecih0(m) = vd(m)
      do l=1,ndim_h0
        vcm(ijm+l) = vu(m,l)
      end do
      ijm = ijm+ndim_h0
    end do
  else
    mmspace = mroot*3+10
    do i=1,mmspace
      indx(i) = (i-1)*ndim_h0
    end do

    do i=1,ndim_h0
      mn = i*(i+1)/2
      vad0(i) = vector2(mn)
    end do

    nxh = ndim_h0*(ndim_h0+1)/2
    nxb = ndim_h0*max_kspace

    call mma_allocate(vb1,max_h0*max_kspace,label='vb1')
    call mma_allocate(vb2,max_h0*max_kspace,label='vb2')
    call basis_2(ndim_h0,vb1,nxb,vad0,vector2,nxh)

    do m=1,mroot
      ecih0(m) = escf(m)
    end do

    kval = mroot*2
    call hymat_2(max_root,max_kspace,ndim_h0,kval,mroot,dcrita,ecih0,vcm,indx,vector2,nxh,vb1,vb2,nxb,vad0)
    !vcm(1:mroot*ndim_h0) = vb1(1:mroot*ndim_h0)
    ! save ci vector in h0 into vb2
    !numh0 = nci_h0 !iw_sta(2,1)
    !vb2(1:numh0*mroot) = Zero
    !if (logic_mr) then
    !  idx1 = 0
    !  idx2 = 0
    !  do i=1,mroot
    !    do j=1,ndim_h0
    !      m = irfno(j)
    !      vb2(idx1+m) = vb1(idx2+j)
    !    end do
    !    idx1 = idx1+numh0
    !    idx2 = idx2+ndim_h0
    !  end do
    !else
    !  vb2(1:numh0*mroot) = vb1(1:numh0*mroot)
    !end if
    call mma_deallocate(vb1)
    call mma_deallocate(vb2)
  end if
  call mma_deallocate(vcm)
end if

write(u6,*)
do m=1,mroot
  write(u6,'(5x,a7,i5,f18.8)') ' root,',m,ecih0(m)
end do
write(u6,*)

!filename = tmpdir(1:len_str)//'/fort7'
!len = len_str+6
!open(nf7,file=filename(1:len),form='unformatted')
!write(nf7) vb2(1:numh0*mroot)
!close(nf7)

!open(100,file='tmp.dat')
!do i=1,ndim_h0
!  write(100,'(1x,f18.9,2i8)') vb1(i),i,irfno(i)
!end do
!write(100,*) 'ndim_h0=',ndim_h0,'v0=',numh0
!idx1 = 0
!idx2 = 0
!do i=1,mroot
!  write(100,*) 'mroot=',mroot
!  do j=1,numh0
!    write(100,'(1x,i2,1x,i5,2x,f18.9)') 1,j,vb2(j)
!    m = ifrno(j)
!    write(100,'(2(1x,i8,1x,f18.9))') j,vb2(j+idx1),m,vb1(m+idx2)
!  end do
!  idx1 = idx1+numh0
!  idx2 = idx2+ndim_h0
!end do
!close(100)
!call abend()

return

end subroutine geth0

subroutine formh0()

use gugaci_global, only: irf, irfno, lenvec, log_prod, logic_mr, max_vector, ndim_h0, vector1, vector2, vint_ci
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp) :: i, iconf1, iconf2, iconfmax, iconfmin, ii, iicc, ir1, ir2, mnh0, mnrh0, num
real(kind=wp), allocatable :: buff(:)

num = ndim_h0*(ndim_h0+1)/2
if (num > max_vector) then
  write(u6,*) ' not enough space to store h0 matrix',num
# ifndef MOLPRO
  call abend()
# endif
  !call abend()
end if

vector2(1:lenvec) = Zero
log_prod = 2
call readint(1,vint_ci)
! act complete loop
call cloop_in_act()
! dbl- act loop
call ploop_in_act()

! save the h0 matrix into sracth file 23 to use later
if (logic_mr) then ! rst
  call mma_allocate(buff,ndim_h0,label='buff')
  buff(1:ndim_h0) = vector1(1:ndim_h0)
  vector1(1:lenvec) = Zero
  mnh0 = ndim_h0*(ndim_h0+1)/2
  vector1(1:mnh0) = vector2(1:mnh0)
  mnrh0 = irf*(irf+1)/2
  vector2(1:mnrh0) = Zero
  do ir1=1,irf
    iconf1 = irfno(ir1)
    ii = ir1*(ir1+1)/2
    vector2(ii) = buff(iconf1)
    do ir2=1,ir1-1
      iconf2 = irfno(ir2)
      ii = ir1*(ir1-1)/2+ir2
      iconfmax = max(iconf1,iconf2)
      iconfmin = min(iconf1,iconf2)
      iicc = iconfmax*(iconfmax-1)/2+iconfmin
      vector2(ii) = vector1(iicc)
    end do
  end do
  ndim_h0 = irf
  call mma_deallocate(buff)
else
  do i=1,ndim_h0   ! rcas
    ii = i*(i-1)/2+i
    vector2(ii) = vector2(ii)+vector1(i)
  end do
end if

!do i=1,ndim_h0
!  write(u6,'(i4,1x,f12.6)') i,vector1(i)
!end do
!call abend()
write(u6,'(a24,i5)') ' dimension of h0 space= ',ndim_h0
log_prod = 1

!open(100,file='h0_new')
!do i=1,num
!  write(100,*) i,vector2(i)
!end do
!close(100)
!call abend()

return

end subroutine formh0

! mroot=2 min_space=4
!    b1=(0,0,1.0,0,0,0...) b2=(0,0,0,0,1.0,0...)

subroutine basis_2(ndim,vb1,nxb,vad,th,nxh)

use gugaci_global, only: ifrno, indx, logic_mr, max_kspace, mjn, mroot
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: ndim, nxb, nxh
real(kind=wp), intent(out) :: vb1(max_kspace*ndim)
real(kind=wp), intent(in) :: vad(ndim), th(nxh)
integer(kind=iwp) :: i, ib, ij, ijh, j, l, m, m0, mief, mjnj
real(kind=wp) :: fenmu, vadi
integer(kind=iwp), allocatable :: ijb1(:)
real(kind=wp), parameter :: dcrita = 1.0e-6_wp, epc = 5.0e-3_wp

vb1(:) = Zero

do j=1,mroot
  ij = indx(j)
  mief = mjn(j)
  if (logic_mr) then
    mjnj = mjn(j)
    mief = ifrno(mjnj)
  end if
  do l=1,ndim
    vb1(ij+l) = Zero
  end do
  vb1(ij+mief) = One
end do

!=======================================================================
call mma_allocate(ijb1,mroot,label='ijb1')
j = mroot
do m=1,mroot
  i = mjn(m)
  if (logic_mr) then
    mjnj = mjn(m)
    i = ifrno(mjnj)
  end if
  vadi = vad(i)
  ijh = i*(i-1)/2
  j = j+1
  ijb1(m) = indx(j)
  do l=1,i-1
    fenmu = vadi-vad(l)
    if (abs(fenmu) < epc) fenmu = epc
    vb1(ijb1(m)+l) = th(ijh+l)/fenmu
  end do
  do l=i+1,ndim
    fenmu = vadi-vad(l)
    if (abs(fenmu) < epc) fenmu = epc
    ijh = l*(l-1)/2
    vb1(ijb1(m)+l) = th(ijh+i)/fenmu
  end do
end do
call mma_deallocate(ijb1)

!-----------------------------------------------------------------------
! write out basis
!-----------------------------------------------------------------------
!write(nf2,*) ' l    vb5       vb6      vb7       vb8'
!
!do l=1,ndim
!  write(nf2,'(2x,i5,4f10.4)') l,vb1(indx(5)+l),vb1(indx(6)+l),vb1(indx(7)+l),vb1(indx(8)+l)
!end do
!call abend()   !wyb_tmp
!-----------------------------------------------------------------------
do m0=1,mroot
  ib = m0+mroot
  call orthnor(ndim,ib,dcrita,vb1,nxb)
end do

return

end subroutine basis_2

subroutine minevalue(iselcsf_occ)

use gugaci_global, only: escf, irf, irfno, logic_mr, LuCiDia, max_innorb, max_ref, mjn, mroot, nci_dim, nci_h0, norb_act, &
                         norb_all, norb_dz, nwalk, vector1, vector2
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(out) :: iselcsf_occ(max_innorb,max_ref)
integer(kind=iwp) :: i, ij, io, j, jm, l, m, ndimh0
integer(kind=iwp), allocatable :: iwalktmp(:)
real(kind=wp) :: am

call read_ml(lucidia,vector1,nci_dim,1)

vector2(1:nci_dim) = vector1(1:nci_dim)
ndimh0 = nci_h0 !iw_sta(2,1)

if (.not. logic_mr) then
  do i=1,mroot
    l = 1
    am = vector2(l)
    do j=1,ndimh0
      if ((vector2(j) /= Zero) .and. (vector2(j) < am)) then
        l = j
        am = vector2(j)
      end if
    end do
    mjn(i) = l
    vector2(l) = Zero
  end do
else
  do i=1,mroot
    l = 1
    am = vector2(l)
    do j=1,irf
      jm = irfno(j)
      if ((vector2(jm) /= Zero) .and. (vector2(jm) < am)) then
        l = jm
        am = vector2(jm)
      end if
    end do
    mjn(i) = l
    vector2(l) = Zero
  end do
end if

do m=1,mroot
  escf(m) = vector1(mjn(m))
end do

call mma_allocate(iwalktmp,norb_all,label='iwalktmp')
write(u6,*) '   mjn(k) :',mroot
do m=1,mroot
  call found_a_config(mjn(m),escf(m),1)
  do i=1,norb_all
    iwalktmp(i) = nwalk(norb_all-i+1)
  end do
  ij = norb_dz
  if (m <= 2*mroot) then
    do io=1,norb_act
      ij = ij+1
      if (iwalktmp(ij) == 3) iselcsf_occ(io,m) = 3
      if (iwalktmp(ij) == 2) iselcsf_occ(io,m) = 2
      if (iwalktmp(ij) == 1) iselcsf_occ(io,m) = 1
      if (iwalktmp(ij) == 0) iselcsf_occ(io,m) = 0
    end do
    !write(u6,'(16(i1))') iwalktmp(norb_dz+1:norb_dz+norb_act)
  end if
  !write(u6 ,'(2x,2i8,f18.8)') m,mjn(m),escf(m)
  !write(nf2,'(2x,2i8,f18.8)') m,mjn(m),escf(m)
end do
write(u6,*)
call mma_deallocate(iwalktmp)

end subroutine minevalue

!subroutine orthnor_ab(n,av,bv,id)  !bv:basis, av:vector for orth a
!
!use Constants, only: Zero
!use Definitions, only: wp, iwp
!
!implicit none
!integer(kind=iwp), intent(in) :: n, id
!real(kind=wp), intent(inout) :: av(n)
!real(kind=wp), intent(in) :: bv(n)
!integer(kind=iwp) :: i
!real(kind=wp) :: s
!real(kind=wp), parameter :: dcrita = 1.0e-10_wp
!real(kind=wp), external :: ddot_
!
!if (id == 0) then
!  ! orthogonalization av,bv
!  s = ddot_(n,av,1,bv,1)
!  do i=1,n
!    av(i) = av(i)-s*bv(i)
!  end do
!end if
!! normalization of av_eigenvector.
!s = Zero
!s = ddot_(n,av,1,av,1)
!s = sqrt(s)
!s = max(s,dcrita)
!do i=1,n
!  av(i) = av(i)/s
!end do
!
!return
!
!end subroutine orthnor_ab

!function ddot_bak(n,dx,dy)
!
!use Constants, only: Zero
!use Definitions, only: wp, iwp
!
!implicit none
!integer(kind=iwp), intent(in) :: n
!real(kind=wp), intent(in) :: dx(n), dy(n)
!integer(kind=iwp) :: l
!real(kind=wp) :: ddot_bak
!real(kind=wp) :: s
!
!s = Zero
!do l=1,n
!  s = s+dx(l)*dy(l)
!end do
!ddot_bak = s
!
!return
!
!end function ddot_bak

!subroutine matrmk_1(k)
!
!use gugaci_global, only: LuCiTv1, LuCiTv2, nci_dim, vector1, vector2, vp
!use Constants, only: Zero
!use Definitions, only: wp, iwp, u6
!
!implicit none
!integer(kind=iwp), intent(in) :: k
!integer(kind=iwp) :: i, ibas, ij, il, jbas, l
!real(kind=wp) :: vsumtmp
!
!do ibas=1,k
!  call read_ml(lucitv1,vector1,nci_dim,ibas)
!  ij = ibas*(ibas-1)/2
!  do jbas=1,ibas
!    call read_ml(lucitv2,vector2,nci_dim,jbas)
!    vsumtmp = Zero
!    do l=1,nci_dim
!      vsumtmp = vsumtmp+vector1(l)*vector2(l)
!    end do
!    vp(ij+jbas) = vsumtmp
!  end do
!end do
!write(u6,*)
!il = 0
!do l=1,k
!  write(u6,1112) (vp(i),i=il+1,il+l)
!  il = il+l
!end do
!write(u6,*)
!
!return
!
!1112 format(2x,20f14.8)
!
!end subroutine matrmk_1
