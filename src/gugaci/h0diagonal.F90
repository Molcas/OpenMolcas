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

subroutine hymat_2(maxroot,minspace,ndim,kval,mroot,dcrita,eval,vcm,indx,th,nxh,vb1,vb2,nxb,vad)
!***********************************************************************
! The sub. based on the algorithm of Davidson [E.R. Davison, J.
! Comput. Phys. 17, 87 (1975)] for searching the m-th eigenvalue and
! eigenvector of the symmetric matrix. This sub. is an improved
! version of the one given by J. Weber, R. Lacroix and G. Wanner
! in Computers & Chemistry, 4, 55 (1980).
!***********************************************************************

use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: maxroot, minspace, ndim, mroot, indx(minspace), nxh, nxb
integer(kind=iwp), intent(inout) :: kval
real(kind=wp), intent(in) :: dcrita, th(nxh), vad(ndim)
real(kind=wp), intent(inout) :: eval(maxroot), vb1(nxb), vb2(nxb)
real(kind=wp), intent(out) :: vcm(ndim*mroot)
integer(kind=iwp) :: i, ijb, ijm, ijmb1, indxm, iterat, jib1, jib2, jicm, k, l, m, mmspace, mn, mrsta, mrsta0, mval, nroot
real(kind=wp) :: depcc, tm, vukm
real(kind=wp), allocatable :: deff(:), ecrita(:), eeval(:), residvb(:), valpha(:), vd(:), ve(:), vp(:), vu(:,:)
real(kind=wp), parameter :: depc = 1.0e-7_wp

!**************************************************************

!write(u6,*) 'generate vector vb2 from matrix a and vector vb1'

!**************************************************************

! debugging
!vb1(1:ndim) = Zero
!vb2(1:ndim) = Zero
!vb1(5) = One
!call abprod2(ndim,1,kval,th,nxh,vb1,vb2,nxb,vad)
!do i=1,ndim
!  write(u6,'(i8,f18.8)') i,vb2(i)
!end do
!call abend()

call mma_allocate(deff,maxroot,label='deff')
call mma_allocate(ecrita,maxroot,label='ecrita')
call mma_allocate(eeval,maxroot,label='eeval')
call mma_allocate(residvb,maxroot,label='residvb')
call mma_allocate(valpha,maxroot,label='valpha')
call mma_allocate(vd,minspace,label='vd')
call mma_allocate(ve,minspace,label='ve')
call mma_allocate(vp,minspace*(minspace+1)/2,label='vp')
call mma_allocate(vu,minspace,minspace,label='vu')

ecrita(:) = 1.0e-8_wp
iterat = 1
mrsta = 1
call abprod2(ndim,1,kval,th,nxh,vb1,vb2,nxb,vad)
call matrmk2(ndim,1,kval,indx,vp,vb1,vb2,nxb)

!=======================================================================
write(u6,*)
l = 0
do k=1,kval
  write(u6,1112) (vp(i),i=l+1,l+k)
  l = l+k
end do
write(u6,*)
!=======================================================================

!write(u6,*) 'diagonalization of matrix p(j,j)'

!=======================================================================
do
  iterat = iterat+1
  if (iterat == 200) then
    write(u6,*) ' h0 space fail to converged'
    write(u6,*) ' program stop'
    call abend()
  end if
  eeval(1:mroot) = eval(1:mroot)
  call hotred(minspace,kval,vp,vd,ve,vu)
  call qlcm(minspace,kval,vd,ve,vu)

  valpha(mrsta:mroot) = Zero
  do m=mrsta,mroot
    eval(m) = vd(m)
    do k=1,kval
      tm = abs(vu(k,m))
      if (valpha(m) < tm) valpha(m) = tm   !max(vu(k,m),k=1,kval)
    end do
    valpha(m) = 1-valpha(m)*valpha(m)
    if (valpha(m) > depc) valpha(m) = sqrt(valpha(m))
    deff(m) = abs(eval(m)-eeval(m))
  end do
  !*********************************************************************
  !
  ! construction of the new approximate eigenvector vcm(n)'
  !
  !*********************************************************************
  ijm = indx(mrsta)
  vcm(ijm+1:ndim*mroot) = Zero
  do k=1,kval
    ijb = indx(k)
    do m=mrsta,mroot
      ijm = indx(m)
      vukm = vu(k,m)
      do i=1,ndim
        vcm(ijm+i) = vcm(ijm+i)+vukm*vb1(ijb+i)
      end do
    end do
  end do

  !do i=1,ndim
  !  write(u6,'(i8,f18.8)') i,vcm(i)
  !end do

  write(u6,1113) iterat,kval,(m,eval(m),deff(m),m=mrsta,mroot)

  nroot = mroot-mrsta+1

  if (kval /= mroot*2) then

    mrsta0 = mrsta
    do m=mrsta0,mroot
      !if ((deff(m) < ecrita(m)) .or. (valpha(m) < dcrita)) then     ! conv
      if ((m == mrsta) .and. (deff(m) < ecrita(m))) then             ! conv
        !if (valpha(m) < dcrita) then                                ! conve
        mrsta = mrsta+1
      end if
    end do
    !mrsta = mrsta+mrooted
    nroot = mroot-mrsta+1
    if (mrsta > mroot) then
      write(u6,*)
      write(u6,*) mroot,' roots are convegenced,after',iterat,' iterat'
      exit
    end if

  end if

  mmspace = min(mroot*3+10,ndim)
  !nroot = 1             ! bbs debug error?
  if (kval+nroot > mmspace) then

    !===== start  reset_basis ==========================================

    do m=mrsta,kval
      indxm = indx(m)
      vb1(indxm+1:indxm+ndim) = Zero
    end do
    do m=mrsta,mroot
      ijm = indx(m)
      do k=1,kval
        ijb = indx(k)
        vukm = vu(k,m)
        do l=1,ndim
          vb1(ijm+l) = vb1(ijm+l)+vukm*vb2(ijb+l)  ! h*cm-->vb1
        end do
      end do
    end do
    ijm = indx(mrsta)
    vb2(ijm+1:ijm+ndim*nroot) = vb1(ijm+1:ijm+ndim*nroot)  ! h*cm-->vb
    vb1(ijm+1:ijm+ndim*nroot) = vcm(ijm+1:ijm+ndim*nroot)  !   cm-->vb
    residvb(mrsta:mroot) = Zero

    mval = mroot
    do m=mrsta,mroot
      ijm = indx(m)
      mval = mval+1
      ijmb1 = indx(mval)
      do l=1,ndim
        depcc = eval(m)-vad(l)
        if (abs(depcc) < depc) depcc = depc
        vb1(ijmb1+l) = (vb2(ijm+l)-vcm(ijm+l)*eval(m))/depcc
        residvb(m) = residvb(m)+vb1(ijmb1+l)*vb1(ijmb1+l)
      end do
      call orthnor(ndim,mval,dcrita,vb1,nxb)
    end do

    vb2(indx(mroot+1)+1:indx(kval+1)) = Zero
    kval = mroot+nroot
    mval = mroot+1
    mn = mroot*(mroot+1)/2
    vp(1:mn) = Zero
    mn = 0
    do m=1,mroot
      mn = mn+m
      vp(mn) = eval(m)
    end do
    !===== end  reset_basis ============================================

  else

    ! form the (j+1)-th approximate vector, vb1(n,j+1)

    jib1 = indx(kval)
    jicm = indx(mrsta)

    residvb(mrsta:mroot) = Zero
    do m=mrsta,mroot
      jib1 = jib1+ndim
      do k=1,kval
        jib2 = indx(k)
        do l=1,ndim
          vb1(jib1+l) = vb1(jib1+l)+vu(k,m)*vb2(jib2+l)
        end do
      end do

      do l=1,ndim
        depcc = eval(m)-vad(l)
        if (abs(depcc) < depc) depcc = depc
        vb1(jib1+l) = (vb1(jib1+l)-vcm(jicm+l)*eval(m))/depcc
        residvb(m) = residvb(m)+vb1(jib1+l)*vb1(jib1+l)
      end do
      jicm = jicm+ndim
    end do
    mval = kval+1
    do m=mrsta,mroot
      kval = kval+1
      call orthnor(ndim,kval,dcrita,vb1,nxb)
    end do

  end if

  call abprod2(ndim,mval,kval,th,nxh,vb1,vb2,nxb,vad)
  call matrmk2(ndim,mval,kval,indx,vp,vb1,vb2,nxb)

  !=====  write out p_matrix ===========================================
  !write(u6,*)
  !write(nf2,*)
  !l = 0
  !do k=1,kval
  !  write(u6,1112) (vp(i),i=l+1,l+k)
  !  write(nf2,1112) (vp(i),i=l+1,l+k)
  !  l = l+k
  !end do
  !write(u6,*)
  !=====  write out p_matrix ===========================================

end do
call mma_deallocate(deff)
call mma_deallocate(ecrita)
call mma_deallocate(eeval)
call mma_deallocate(residvb)
call mma_deallocate(valpha)
call mma_deallocate(vd)
call mma_deallocate(ve)
call mma_deallocate(vp)
call mma_deallocate(vu)

! copy ci vector to VB1
do m=1,mroot
  ijm = indx(m)
  vb1(ijm+1:ijm+ndim) = vcm(ijm+1:ijm+ndim)
  !write(nf1) eval(m),(vcm(ijm+i),i=1,ndim)
  !write(u6,'(5(1x,f18.9),1x,i2)') (vcm(ijm+i),i=1,4),vcm(35)
end do

return

1112 format(2x,20f14.8)
1113 format(2i3,10(2x,i2,f14.8,f14.8))

end subroutine hymat_2

subroutine matrmk2(n,k1,k2,indx,p,vb1,vb2,nxb)

use Constants, only: Zero
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: n, k1, k2, indx(30), nxb
real(kind=wp), intent(inout) :: p(465)
real(kind=wp), intent(in) :: vb1(nxb), vb2(nxb)
integer(kind=iwp) :: i, iijj, ij, j, ji, l

do i=k1,k2
  iijj = i*(i-1)/2
  ij = indx(i)
  do j=1,i
    ji = indx(j)
    p(iijj+j) = Zero
    !-------------------------------------------------------------------
    do l=1,n
      !-----------------------------------------------------------------
      p(iijj+j) = p(iijj+j)+vb1(ij+l)*vb2(ji+l)
    end do
  end do
end do

return

end subroutine matrmk2

subroutine abprod2(n,k1,k2,th,nxh,vb1,vb2,nxb,vad)

use gugaci_global, only: indx
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: n, k1, k2, nxh, nxb
real(kind=wp), intent(in) :: th(nxh), vb1(nxb), vad(n)
real(kind=wp), intent(inout) :: vb2(nxb)
integer(kind=iwp) :: i, ij, j, l, mn

ij = 0
do j=k1,k2
  ij = indx(j)
  do i=1,n
    vb2(ij+i) = vad(i)*vb1(ij+i)
  end do
end do
!-----------------------------------------------------------------------
do i=2,n
  mn = i*(i-1)/2
  do j=k1,k2
    ij = indx(j)
    do l=1,i-1
      vb2(ij+i) = vb2(ij+i)+th(mn+l)*vb1(ij+l)
      vb2(ij+l) = vb2(ij+l)+th(mn+l)*vb1(ij+i)
    end do
  end do
end do
!-----------------------------------------------------------------------
return

end subroutine abprod2

subroutine orthnor(n,j,dcrita,vb1,nxb)

use gugaci_global, only: indx
use Constants, only: Zero
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: n, j, nxb
real(kind=wp), intent(in) :: dcrita
real(kind=wp), intent(inout) :: vb1(nxb)
integer(kind=iwp) :: i, ij, ji, jm, l
real(kind=wp) :: s, smax1, smax2

ji = indx(j)
if (j /= 1) then
  jm = j-1
  smax2 = huge(smax2)
  do
    smax1 = Zero
    do l=1,jm
      s = Zero
      ij = indx(l)
      do i=1,n
        s = s+vb1(ij+i)*vb1(ji+i)
      end do
      smax1 = max(smax1,abs(s))
      do i=1,n
        vb1(ji+i) = vb1(ji+i)-s*vb1(ij+i)
      end do
    end do

    if (smax1 < dcrita) exit
    if (smax1 > smax2) then
      write(u6,*) 'diagonalization procedure is non-convergent.'
#     ifndef MOLPRO
      call abend()
#     endif
      !call abend()
    end if
    smax2 = smax1
  end do
end if
! normalization of j-th eigenvector.
s = Zero
do i=1,n
  s = s+vb1(ji+i)*vb1(ji+i)
end do
s = sqrt(s)
do i=1,n
  vb1(ji+i) = vb1(ji+i)/s
end do

return

end subroutine orthnor

subroutine norm_a(n,av)  !bv:basis, av:vector for orth and norm

use Constants, only: Zero
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: n
real(kind=wp), intent(inout) :: av(n)
integer(kind=iwp) :: i
real(kind=wp) :: s
real(kind=wp), parameter :: dcrita = 1.0e-10_wp
real(kind=wp), external :: ddot_

! normalization of av_eigenvector.
s = Zero
s = ddot_(n,av,1,av,1)
s = sqrt(s)
s = max(s,dcrita)
do i=1,n
  av(i) = av(i)/s
end do

return

end subroutine norm_a

subroutine orth_ab(n,av,bv)  !bv:basis, av:vector for orth

use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: n
real(kind=wp), intent(inout) :: av(n)
real(kind=wp), intent(in) :: bv(n)
integer(kind=iwp) :: i
real(kind=wp) :: s
real(kind=wp), external :: ddot_

! orthogonalization av,bv
s = ddot_(n,av,1,bv,1)

do i=1,n
  av(i) = av(i)-s*bv(i)
end do

return

end subroutine orth_ab
