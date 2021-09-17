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

implicit none
!integer, parameter :: maxroot = 10, minspace = 40)
integer :: maxroot, minspace, ndim, kval, mroot, indx(minspace), nxh, nxb
real*8 :: dcrita, eval(maxroot), vcm(ndim*mroot), th(nxh), vb1(nxb), vb2(nxb), vad(ndim)
integer :: i, ijb, ijm, ijmb1, indxm, iterat, jib1, jib2, jicm, k, l, m, mmspace, mn, mrsta, mrsta0, mval, nroot
real*8 :: deff(maxroot), depcc, ecrita(maxroot), eeval(maxroot), residvb(maxroot), tm, valpha(maxroot), vd(minspace), &
          ve(minspace), vp(minspace*(minspace+1)/2), vu(minspace,minspace), vukm
real*8, parameter :: depc = 1.0d-7

!**************************************************************

!write(6,*) 'generate vector vb2 from matrix a and vector vb1'

!**************************************************************

! debugging
!vb1(1:ndim) = 0.d0
!vb2(1:ndim) = 0.d0
!vb1(5) = 1.d0
!call abprod2(ndim,1,kval,th,nxh,vb1,vb2,nxb,vad)
!do i=1,ndim
!  write(6,'(i8,f18.8)') i,vb2(i)
!end do
!call abend()

ecrita = 1.e-8
iterat = 1
mrsta = 1
call abprod2(ndim,1,kval,th,nxh,vb1,vb2,nxb,vad)
call matrmk2(ndim,1,kval,indx,vp,vb1,vb2,nxb)

!=======================================================================
write(6,*)
l = 0
do k=1,kval
  write(6,1112) (vp(i),i=l+1,l+k)
  l = l+k
end do
write(6,*)
!=======================================================================

!write(6,*) 'diagonalization of matrix p(j,j)'

!=======================================================================
do
  iterat = iterat+1
  if (iterat == 200) then
    write(6,*) ' h0 space fail to converged'
    write(6,*) ' program stop'
    call abend()
  end if
  eeval(1:mroot) = eval(1:mroot)
  call hotred(minspace,kval,vp,vd,ve,vu)
  call qlcm(minspace,kval,vd,ve,vu)

  valpha(mrsta:mroot) = 0.d0
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
  vcm(ijm+1:ndim*mroot) = 0.0d0
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
  !  write(6,'(i8,f18.8)') i,vcm(i)
  !end do

  write(6,1113) iterat,kval,(m,eval(m),deff(m),m=mrsta,mroot)
  1113 format(2i3,10(2x,i2,f14.8,f14.8))

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
      write(6,*)
      write(6,*) mroot,' roots are convegenced,after',iterat,' iterat'
      exit
    end if

  end if

  mmspace = min(mroot*3+10,ndim)
  !nroot = 1             ! bbs debug error?
  if (kval+nroot > mmspace) then

    !===== start  reset_basis ==========================================

    do m=mrsta,kval
      indxm = indx(m)
      vb1(indxm+1:indxm+ndim) = 0.0d0
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
    residvb(mrsta:mroot) = 0.d0

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

    vb2(indx(mroot+1)+1:indx(kval+1)) = 0.d0
    kval = mroot+nroot
    mval = mroot+1
    mn = mroot*(mroot+1)/2
    vp(1:mn) = 0.0d0
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

    residvb(mrsta:mroot) = 0.d0
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
  !write(6,*)
  !write(nf2,*)
  !l = 0
  !do k=1,kval
  !  write(6,1112) (vp(i),i=l+1,l+k)
  !  write(nf2,1112) (vp(i),i=l+1,l+k)
  !  l = l+k
  !end do
  !write(6,*)
  !=====  write out p_matrix ===========================================

end do

! copy ci vector to VB1
do m=1,mroot
  ijm = indx(m)
  vb1(ijm+1:ijm+ndim) = vcm(ijm+1:ijm+ndim)
  !write(nf1) eval(m),(vcm(ijm+i),i=1,ndim)
  !write(6,'(5(1x,f18.9),1x,i2)') (vcm(ijm+i),i=1,4),vcm(35)
end do

return

1112 format(2x,20f14.8)

end subroutine hymat_2

subroutine matrmk2(n,k1,k2,indx,p,vb1,vb2,nxb)

implicit none
integer :: n, k1, k2, indx(30), nxb
real*8 :: p(465), vb1(nxb), vb2(nxb)
integer :: i, iijj, ij, j, ji, l

do i=k1,k2
  iijj = i*(i-1)/2
  ij = indx(i)
  do j=1,i
    ji = indx(j)
    p(iijj+j) = 0.0d0
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

implicit none
integer :: n, k1, k2, nxh, nxb
real*8 :: th(nxh), vb1(nxb), vb2(nxb), vad(n)
integer :: i, ij, j, l, mn
!real*8, allocatable :: buff(:)

!allocate(buff(n))
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

implicit none
integer :: n, j, nxb
real*8 :: dcrita, vb1(nxb)
integer :: i, ij, ji, jm, l
real*8 :: s, smax1, smax2

ji = indx(j)
if (j /= 1) then
  jm = j-1
  smax2 = 1.d10
  do
    smax1 = 0.0d0
    do l=1,jm
      s = 0.0d0
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
      write(6,*) 'diagonalization procedure is non-convergent.'
#     ifndef MOLPRO
      call abend()
#     endif
      !call abend()
    end if
    smax2 = smax1
  end do
end if
! normalization of j-th eigenvector.
s = 0.0d0
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

implicit none
integer :: n
real*8 :: av(n)
integer :: i
real*8 :: s
real*8, parameter :: dcrita = 1.0d-10
real*8, external :: ddot_

! normalization of av_eigenvector.
s = 0.0d0
s = ddot_(n,av,1,av,1)
s = sqrt(s)
s = max(s,dcrita)
do i=1,n
  av(i) = av(i)/s
end do

return

end subroutine norm_a

subroutine orth_ab(n,av,bv)  !bv:basis, av:vector for orth

implicit none
integer :: n
real*8 :: av(n), bv(n)
integer :: i
real*8 :: s
real*8, external :: ddot_

! orthogonalization av,bv
s = ddot_(n,av,1,bv,1)

do i=1,n
  av(i) = av(i)-s*bv(i)
end do

return

end subroutine orth_ab
