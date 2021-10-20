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

!module_grad works for the ci gradient calculations.
! subroutine ci_grad
! subroutine ci_density_label
! subroutine convert_vector
! subroutine trans_ijkl_intpos
! subroutine trans_intpos_ijkl
! subroutine square_canonical
! subroutine moread
! subroutine lagran
! subroutine cidensity_test

!subroutine ci_grad()
!
!implicit none
!
!!================================================
!! the main subroutine for ci gradient calculations.
!! ican_a and ican_b save canonical order.
!! e(*) save the scf orbital energies.
!! xlgrn(*) save the lagrangian matrix.
!
!! ci gradient is still in development, so we do not include it in current
!! version of xian-ci
!
!return
!
!!write(u6,'(a18,2x,f10.2,2x,a1)') 'end of grad, takes',sc5-sc0,'s'
!end subroutine ci_grad

!subroutine convert_vector()
!
!implicit none
!
!!=====================================================
!! this just uses at debug.
!
!!open(10,file='cigmsvector')
!!do i=1,nci_dim
!!  read(10,*) iiii,val1,j,val2
!!  vector1(j) = val1
!!  if ((val2 > Zero) .and. (val1 < Zero)) vector1(j) = -vector1(j)
!!  if ((val2 < Zero) .and. (val1 > Zero)) vector1(j) = -vector1(j)
!!end do
!!close(10)
!
!end subroutine convert_vector

!subroutine moread(ii,jj,kk,ll,val)
!
!use gugaci_global, only: ican_a, ican_b, vector1
!use Definitions, only: wp, iwp
!
!implicit none
!integer(kind=iwp), intent(in) :: ii, jj, kk, ll
!real(kind=wp), intent(in) :: val
!integer(kind=iwp) :: lri, lrj, lrk, lrl, lrn, nij, nijkl, nkl
!
!!====================================================
!! transfer the ci mo indices to scf mo indices
!! and use library function ishft() to save
!! two electron mo indices to dm2_index().
!
!lrj = min(ii,jj)
!lri = max(ii,jj)
!lrl = min(kk,ll)
!lrk = max(kk,ll)
!if (lri < lrk) then
!  lrn = lrk
!  lrk = lri
!  lri = lrn
!  lrn = lrl
!  lrl = lrj
!  lrj = lrn
!end if
!if ((lri == lrk) .and. (lrj < lrl)) then
!  lrn = lrj
!  lrj = lrl
!  lrl = lrn
!end if
!
!nij = ican_a(lri)+lrj
!nkl = ican_a(lrk)+lrl
!nijkl = ican_b(nij)+nkl
!vector1(nijkl) = val
!!write(nf2,'(4i4,f18.10)') lri,lrj,lrk,lrl,vector1(nijkl)
!
!end subroutine moread

subroutine trans_ijkl_intpos(ii,jj,kk,ll,nxo)

use gugaci_global, only: ican_a, ican_b
use Definitions, only: iwp

implicit none
integer(kind=iwp), intent(in) :: ii, jj, kk, ll
integer(kind=iwp), intent(out) :: nxo
integer(kind=iwp) :: lri, lrj, lrk, lrl, lrn, nij, nkl

!====================================================
! transfer the ci mo indices to scf mo indices
! and use library function ishft() to save
! two electron mo indices to dm2_index().
!write(nf2,'(4i4)') ii,jj,kk,ll

!i = map_orb_order_t(ii)
!j = map_orb_order_t(jj)
!k = map_orb_order_t(kk)
!l = map_orb_order_t(ll)

lrj = min(ii,jj)
lri = max(ii,jj)
lrl = min(kk,ll)
lrk = max(kk,ll)
if (lri < lrk) then
  lrn = lrk
  lrk = lri
  lri = lrn
  lrn = lrl
  lrl = lrj
  lrj = lrn
end if
if ((lri == lrk) .and. (lrj < lrl)) then
  lrn = lrj
  lrj = lrl
  lrl = lrn
end if

nij = ican_a(lri)+lrj
nkl = ican_a(lrk)+lrl
nxo = ican_b(nij)+nkl

!if (nxo == 3003) then
!  write(nf2,'(8i4,i8)') ii,jj,kk,ll,lri,lrj,lrk,lrl,nxo
!end if

end subroutine trans_ijkl_intpos

!subroutine trans_intpos_ijkl()
!
!implicit none
!
!!=========================================================
!! read the saved two electron mo indices.
!
!return
!
!end subroutine trans_intpos_ijkl

!subroutine lagran_act(x1e)
!
!use gugaci_global, only: dm1, ican_a, ican_b, norb_all, vector1, vector2, xlgrn
!use stdalloc, only: mma_allocate, mma_deallocate
!use Constants, only: Zero, Two
!use Definitions, only: wp, iwp
!
!implicit none
!real(kind=wp), intent(in) :: x1e(50000)
!integer(kind=iwp) :: i, i0, j, j0, k, k0, kl, l, l0, m, mik, mjk, nik, nil, niljk, nimkl, nji, njikl, njk, njmkl, nkl, norbf
!real(kind=wp) :: dum, dumtmp
!real(kind=wp), allocatable :: fock(:,:)
!
!!================================================
!! lyb
!! xlgrn(norb_all,norb_all) is the lagrange matrix.
!
!norbf = n_frz+1
!!norbf = 1
!
!xlgrn(1:n_all,1:n_all) = Zero
!
!call mma_allocate(fock,n_all,n_all,label='fock')
!call lagran_fock(x1e,fock)
!
!! form two electron contributions to the lagrangian with frozen mo
!
!do i=ndbl+1,n_all
!  do j=1,n_frz
!    dum = fock(j,i)*Two
!    do k=norbf,n_all
!      do l=norbf,n_all
!        i0 = max(i,j)
!        j0 = min(i,j)
!        nji = ican_a(i0)+j0
!        k0 = max(k,l)
!        l0 = min(k,l)
!        nkl = ican_a(k0)+l0
!        if (nji >= nkl) then
!          njikl = ican_b(nji)+nkl
!        else
!          njikl = ican_b(nkl)+nji
!        end if
!        dum = dum+dm1(k0,l0)*vector1(njikl)*Two
!
!        i0 = max(i,l)
!        l0 = min(i,l)
!        nil = ican_a(i0)+l0
!        j0 = max(j,k)
!        k0 = min(j,k)
!        njk = ican_a(j0)+k0
!        if (nil >= njk) then
!          niljk = ican_b(nil)+njk
!        else
!          niljk = ican_b(njk)+nil
!        end if
!        dum = dum-dm1(k0,l0)*vector1(niljk)
!      end do
!    end do
!    xlgrn(i,j) = xlgrn(i,j)+dum
!  end do
!end do
!
!do i=1,n_frz
!  do j=ndbl+1,n_all
!    dum = fock(j,i)*Two
!    do k=norbf,n_all
!      do l=norbf,n_all
!        i0 = max(i,j)
!        j0 = min(i,j)
!        nji = ican_a(i0)+j0
!        k0 = max(k,l)
!        l0 = min(k,l)
!        nkl = ican_a(k0)+l0
!        if (nji >= nkl) then
!          njikl = ican_b(nji)+nkl
!        else
!          njikl = ican_b(nkl)+nji
!        end if
!        dum = dum+dm1(k0,l0)*vector1(njikl)*Two
!
!        i0 = max(i,l)
!        l0 = min(i,l)
!        nil = ican_a(i0)+l0
!        j0 = max(j,k)
!        k0 = min(j,k)
!        njk = ican_a(j0)+k0
!        if (nil >= njk) then
!          niljk = ican_b(nil)+njk
!        else
!          niljk = ican_b(njk)+nil
!        end if
!        dum = dum-dm1(k0,l0)*vector1(niljk)
!      end do
!    end do
!    xlgrn(i,j) = xlgrn(i,j)+dum
!  end do
!end do
!
!! form two electron contributions to the lagrangian with active mo
!
!do i=norbf,ndbl
!  do j=ndbl+1,n_all
!    dum = Zero
!    do m=norbf,n_all
!      i0 = ican_a(max(i,m))+min(i,m)
!      j0 = ican_a(max(j,m))+min(j,m)
!      kl = 0
!      do k=norbf,n_all
!        dumtmp = Zero
!        do l=norbf,k-1
!          kl = kl+1
!          if (kl > i0) then
!            nimkl = ican_b(kl)+i0
!          else
!            nimkl = ican_b(i0)+kl
!          end if
!
!          if (kl > j0) then
!            njmkl = ican_b(kl)+j0
!          else
!            njmkl = ican_b(j0)+kl
!          end if
!          dumtmp = dumtmp+vector1(nimkl)*vector2(njmkl)
!          !if ((i == 1) .and. (j == 5)) write(nf2,'(6i4,i8,f18.10)') i,j,m,k,l,kl,njmkl,vector2(njmkl)
!        end do
!
!        dum = dum+dumtmp*Two
!
!        kl = kl+1
!        if (kl > i0) then
!          nimkl = ican_b(kl)+i0
!        else
!          nimkl = ican_b(i0)+kl
!        end if
!
!        if (kl > j0) then
!          njmkl = ican_b(kl)+j0
!        else
!          njmkl = ican_b(j0)+kl
!        end if
!
!        dum = dum+vector1(nimkl)*vector2(njmkl)
!        !if (nimkl /= njmkl) write(nf2,'(2i8,2f18.10)') nimkl,njmkl,vector1(nimkl),vector2(njmkl)
!        !if ((i == 1) .and. (j == 5)) write(nf2,'(6i4,i8,f18.10)') i,j,m,k,l,kl,njmkl,vector2(njmkl
!      end do
!    end do
!
!    xlgrn(i,j) = xlgrn(i,j)+dum*Two
!    !write(2,'(a7,2i4,f18.10)') 'xlgrn_2',i,j,xlgrn(i,j)
!  end do
!end do
!do i=norbf,ndbl
!  do j=ndbl+1,n_all
!    dum = Zero
!    do k=norbf,norb_all
!      mik = max(i,k)
!      nik = min(i,k)
!      mjk = max(j,k)
!      njk = min(j,k)
!      dum = dum+fock(mik,nik)*dm1(mjk,njk)
!
!      !write(nf2,'(2f18.10)') x1e(mnik),dm1(mjk,njk)
!    end do
!    xlgrn(i,j) = xlgrn(i,j)+dum
!    !write(2,'(a9,2i4,f18.10)') 'xlgrn_all',i,j,xlgrn(i,j)
!  end do
!end do
!
!do i=ndbl+1,n_all
!  do j=norbf,ndbl
!    dum = Zero
!    do m=norbf,n_all
!      i0 = ican_a(max(i,m))+min(i,m)
!      j0 = ican_a(max(j,m))+min(j,m)
!      kl = 0
!      do k=norbf,n_all
!        dumtmp = Zero
!        do l=norbf,k-1
!          kl = kl+1
!          if (kl > i0) then
!            nimkl = ican_b(kl)+i0
!          else
!            nimkl = ican_b(i0)+kl
!          end if
!
!          if (kl > j0) then
!            njmkl = ican_b(kl)+j0
!          else
!            njmkl = ican_b(j0)+kl
!          end if
!          dumtmp = dumtmp+vector1(nimkl)*vector2(njmkl)
!          !if ((i == 1) .and. (j == 5)) write(nf2,'(6i4,i8,f18.10)') i,j,m,k,l,kl,njmkl,vector2(njmkl)
!        end do
!
!        dum = dum+dumtmp*Two
!
!        kl = kl+1
!        if (kl > i0) then
!          nimkl = ican_b(kl)+i0
!        else
!          nimkl = ican_b(i0)+kl
!        end if
!
!        if (kl > j0) then
!          njmkl = ican_b(kl)+j0
!        else
!          njmkl = ican_b(j0)+kl
!        end if
!
!        dum = dum+vector1(nimkl)*vector2(njmkl)
!        !if (nimkl /= njmkl) write(nf2,'(2i8,2f18.10)') nimkl,njmkl,vector1(nimkl),vector2(njmkl)
!        !if ((i == 1) .and. (j == 5)) write(nf2,'(6i4,i8,f18.10)') i,j,m,k,l,kl,njmkl,vector2(njmkl
!      end do
!    end do
!
!    xlgrn(i,j) = xlgrn(i,j)+dum*Two
!    !write(2,'(a7,2i4,f18.10)') 'xlgrn_2',i,j,xlgrn(i,j)
!  end do
!end do
!
!do i=ndbl+1,n_all
!  do j=norbf,ndbl
!    dum = Zero
!    do k=norbf,norb_all
!      mik = max(i,k)
!      nik = min(i,k)
!      mjk = max(j,k)
!      njk = min(j,k)
!      dum = dum+fock(mik,nik)*dm1(mjk,njk)
!
!      !write(nf2,'(2f18.10)') x1e(mnik),dm1(mjk,njk)
!    end do
!    xlgrn(i,j) = xlgrn(i,j)+dum
!    !write(2,'(a9,2i4,f18.10)') 'xlgrn_all',i,j,xlgrn(i,j)
!  end do
!end do
!call mma_deallocate(fock)
!
!end subroutine lagran_act

!subroutine lagran_all(x1e)
!
!use gugaci_global, only: dm1, ican_a, ican_b, norb_all, vector1, vector2, xlgrn
!use Constants, only: Zero, Two
!use Definitions, only: wp, iwp
!
!implicit none
!real(kind=wp), intent(in) :: x1e(50000)
!integer(kind=iwp) :: i, i0, j, j0, k, kl, l, m, mik, mjk, mnik, nik, nimkl, njk, njmkl
!real(kind=wp) :: dum, dumtmp
!
!!================================================
!! lyb
!! xlgrn(norb_all,norb_all) is the lagrange matrix.
!
!!norbf = n_frz+1
!!norbf = 1
!
!xlgrn(1:norb_all,1:norb_all) = Zero
!
!! form two electron contributions to the lagrangian
!do i=1,norb_all
!  do j=1,norb_all
!    dum = Zero
!    do m=1,norb_all
!      i0 = ican_a(max(i,m))+min(i,m)
!      j0 = ican_a(max(j,m))+min(j,m)
!      kl = 0
!      do k=1,norb_all
!        dumtmp = Zero
!        do l=1,k-1
!          kl = kl+1
!          if (kl > i0) then
!            nimkl = ican_b(kl)+i0
!          else
!            nimkl = ican_b(i0)+kl
!          end if
!
!          if (kl > j0) then
!            njmkl = ican_b(kl)+j0
!          else
!            njmkl = ican_b(j0)+kl
!          end if
!          dumtmp = dumtmp+vector1(nimkl)*vector2(njmkl)
!          !if ((i == 1) .and. (j == 5)) write(nf2,'(6i4,i8,f18.10)') i,j,m,k,l,kl,njmkl,vector2(njmkl)
!        end do
!
!        dum = dum+dumtmp*Two
!
!        kl = kl+1
!        if (kl > i0) then
!          nimkl = ican_b(kl)+i0
!        else
!          nimkl = ican_b(i0)+kl
!        end if
!
!        if (kl > j0) then
!          njmkl = ican_b(kl)+j0
!        else
!          njmkl = ican_b(j0)+kl
!        end if
!
!        dum = dum+vector1(nimkl)*vector2(njmkl)
!        !if (nimkl /= njmkl) write(nf2,'(2i8,2f18.10)') nimkl,njmkl,vector1(nimkl),vector2(njmkl)
!        !if ((i == 1) .and. (j == 5)) write(nf2,'(6i4,i8,f18.10)') i,j,m,k,l,kl,njmkl,vector2(njmkl)
!      end do
!    end do
!
!    xlgrn(i,j) = xlgrn(i,j)+dum*Two
!    !write(2,'(a7,2i4,f18.10)') 'xlgrn_2',i,j,xlgrn(i,j)
!  end do
!end do
!do i=1,norb_all
!  do j=1,norb_all
!    dum = Zero
!    do k=1,norb_all
!      mik = max(i,k)
!      nik = min(i,k)
!      mjk = max(j,k)
!      njk = min(j,k)
!      mnik = ican_a(mik)+nik
!      dum = dum+x1e(mnik)*dm1(mjk,njk)
!
!      !write(nf2,'(2f18.10)') x1e(mnik),dm1(mjk,njk)
!    end do
!    xlgrn(i,j) = xlgrn(i,j)+dum
!    !write(2,'(a9,2i4,f18.10)') 'xlgrn_all',i,j,xlgrn(i,j)
!  end do
!end do
!
!end subroutine lagran_all

!subroutine writedm2(nx)
!
!!use gugaci_global, only: len_str, tmpdir
!use Definitions, only: iwp
!
!implicit none
!integer(kind=iwp), intent(in) :: nx
!!character(len=256) :: filename
!
!!filename = tmpdir(1:len_str)//'/density'
!!len = len_str+8
!! Need debug
!!open(nf22,file=filename(1:len),form='unformatted')
!
!!open(nf22,file='density',form='unformatted')
!!write(nf22) (vector2(i),i=1,nx)
!!close(nf22)
!
!end subroutine writedm2
!
!subroutine readdm2(nx)
!
!!use gugaci_global, only: len_str, tmpdir
!use gugaci_global, only: vector2
!use Constants, only: Zero
!use Definitions, only: iwp
!
!implicit none
!integer(kind=iwp), intent(in) :: nx
!!character(len=256) :: filename
!
!vector2(1:nx) = Zero
!
!!filename = tmpdir(1:len_str)//'/density'
!!len = len_str+8
!
!!open(nf22,file=filename(1:len),form='unformatted')
!
!!open(nf22,file='density',form='unformatted')
!!read(nf22) (vector2(i),i=1,nx)
!!close(nf22)
!
!end subroutine readdm2

!subroutine backtrans_test()
!
!use gugaci_global, only: cf, ican_a, ican_b, norb_all, norb_frz, vector2
!use Constants, only: Zero
!use Definitions, only: wp, iwp
!
!implicit none
!integer(kind=iwp) :: i, i0, j, j0, k, k0, l, l0, nij, nijkl, nkl
!real(kind=wp) :: sum_, val
!
!!-------------------------------------------------------------
!! this subroutine for test the backtrans result by just one density matr
!! ao, such as (ij|kl)
!
!sum_ = Zero
!i0 = 6
!j0 = 6
!k0 = 6
!l0 = 6
!
!do i=norb_frz+1,norb_all
!  do j=norb_frz+1,norb_all
!    if (i >= j) then
!      nij = ican_a(i)+j
!    else
!      nij = ican_a(j)+i
!    end if
!
!    do k=norb_frz+1,norb_all
!      do l=norb_frz+1,norb_all
!        if (k >= l) then
!          nkl = ican_a(k)+l
!        else
!          nkl = ican_a(l)+k
!        end if
!        if (nij >= nkl) then
!          nijkl = ican_b(nij)+nkl
!        else
!          nijkl = ican_b(nkl)+nij
!        end if
!        val = vector2(nijkl)
!        sum_ = sum_+val*cf(i0,i)*cf(j0,j)*cf(k0,k)*cf(l0,l)
!      end do
!    end do
!  end do
!end do
!
!!write(nf2,*) 'num 1 sum=',sum_
!
!end subroutine backtrans_test

!subroutine grad_two()
!
!use gugaci_global, only: dxyz, ican_a, ican_b, numat, vector1 !, len_str, tmpdir
!use stdalloc, only: mma_allocate, mma_deallocate
!use Constants, only: Zero, One, Two, auTokcalmol
!use Definitions, only: wp, iwp, u6
!
!implicit none
!integer(kind=iwp) :: i, i0, ind, j, j0, k, k0, l, l0, ncon, nij, nijkl, nkl, nnij, npat
!real(kind=wp) :: aa, bb, val, val1, val2
!!character(len=256) :: filename
!integer(kind=iwp), allocatable :: index_atom(:,:), ndi0(:), ndj0(:), ndk0(:), ndl0(:)
!real(kind=wp), allocatable :: daoint1(:), daoxyz(:,:), dgxyz(:,:)
!
!npat = numat*(numat+1)/2
!call mma_allocate(index_atom,3,npat,label='index_atom')
!index_atom(1:3,1:npat) = 0
!
!call mma_allocate(ndi0,ndao,label='ndi0')
!call mma_allocate(ndj0,ndao,label='ndj0')
!call mma_allocate(ndk0,ndao,label='ndk0')
!call mma_allocate(ndl0,ndao,label='ndl0')
!call mma_allocate(daoint1,ndao,label='daoint1')
!ndi0(:) = 0
!ndj0(:) = 0
!ndk0(:) = 0
!ndl0(:) = 0
!daoint1(:) = Zero
!
!!filename=  tmpdir(1:len_str)//'/daoints'
!!len = len_str+8
!
!!open(40,file=filename(1:len),form='formatted')
!
!!open(40,file='daoints',form='formatted')
!!read(40,*)
!
!!do i=1,3
!!  read(40,*) (index_atom(i,j),j=1,npat)
!!end do
!!read(40,*) (ndi0(i),i=1,ndao)
!!read(40,*) (ndj0(i),i=1,ndao)
!!read(40,*) (ndk0(i),i=1,ndao)
!!read(40,*) (ndl0(i),i=1,ndao)
!!read(40,*) (daoint1(i),i=1,ndao)
!!
!!close(40)
!
!ncon=0
!
!call mma_allocate(daoxyz,3,numat,label='daoxyz')
!do i=1,numat
!  do j=1,i-1
!    nnij = ican_a(i)+j
!    !write(nf2,*) i,ican_a(i)
!    do k=1,3
!      ind = index_atom(k,nnij)
!      val = Zero
!      do l=1,ind
!        ncon = ncon+1
!        i0 = ndi0(ncon)
!        j0 = ndj0(ncon)
!        k0 = ndk0(ncon)
!        l0 = ndl0(ncon)
!        val2 = daoint1(ncon)
!        !write(2,'(4i4,f18.10)') i0,j0,k0,l0,val2
!        nij = ican_a(j0)+i0
!        nkl = ican_a(l0)+k0
!        nijkl = ican_b(nij)+nkl
!        val1 = vector1(nijkl)
!        aa = One
!        bb = One
!        if (i0 /= j0) aa = Two*aa
!        if (k0 /= l0) bb = Two*bb
!        !===============================================================
!        ! this place should multiple two, because the gradient
!        ! integral always exist the relation that :
!        ! index nij=i*(i-1)/2+j
!        ! index nkl=k*(k-1)/2+l
!        ! but here always the nij >= nkl.
!
!        val = val+val1*val2*aa*bb*Two
!      end do
!      daoxyz(k,i) = daoxyz(k,i)+val
!      daoxyz(k,j) = daoxyz(k,j)-val
!
!      dxyz(k,i) = dxyz(k,i)+val
!      dxyz(k,j) = dxyz(k,j)-val
!    end do
!  end do
!end do
!call mma_deallocate(index_atom)
!call mma_deallocate(ndi0)
!call mma_deallocate(ndj0)
!call mma_deallocate(ndk0)
!call mma_deallocate(ndl0)
!call mma_deallocate(daoint1)
!
!call mma_allocate(dgxyz,3,numat,label='dgxyz')
!!dgxyz(1:3,1:numat) = Zero
!!do i=1,numat
!!  do j=1,3
!!    dgxyz(j,i) = daoxyz(j,i)*auTokcalmol
!!  end do
!!end do
!call mma_deallocate(daoxyz)
!
!
!!write(nf2,'(//10x,'cartesian coordinate derivatives',//3x,'number  atom ',5x,'x',12x,'y',12x,'z',/)')
!
!!do i=1,numat
!!  write(nf2,'(6x,i6,3f13.6)') i,(dgxyz(j,i),j=1,3)
!!end do
!
!dgxyz(1:3,1:numat) = Zero
!do i=1,numat
!  do j=1,3
!    dgxyz(j,i) = dxyz(j,i)*auTokcalmol
!  end do
!end do
!
!!write(nf2,'(//10x,'cartesian coordinate derivatives',//3x,'number  atom ',5x,'x',12x,'y',12x,'z',/)')
!
!do i=1,numat
!  write(u6,'(6x,i6,3f13.6)') i,(dgxyz(j,i),j=1,3)
!end do
!call mma_deallocate(dgxyz)
!
!end subroutine grad_two

!subroutine grad_one_ao()
!
!use gugaci_global, only: dxyz, ican_a, naorbs, numat, p !, len_str, tmpdir
!use stdalloc, only: mma_allocate, mma_deallocate
!use Constants, only: Zero, Two, auTokcalmol
!use Definitions, only: wp, iwp, u6
!
!implicit none
!integer(kind=iwp) :: i, j, k, l, nkl
!!character(len=256) :: filename
!real(kind=wp), allocatable :: dgxyz(:,:), dm1_act(:,:), dmo1xyz(:,:), dsaos(:,:,:)
!
!call mma_allocate(dsaos,3,numat,naorbs*(naorbs+1)/2,label='dsaos')
!dsaos(1:3,1:numat,1:naorbs*(naorbs+1)/2) = Zero
!
!!filename = tmpdir(1:len_str)//'/dfock1'
!!len = len_str+7
!
!!open(500,file=filename(1:len),form='unformatted')
!
!!open(500,file='dfock1',form='unformatted')
!!do i=1,numat
!!  do k=1,3
!!    read(500) (dsaos(k,i,j),j=1,naorbs*(naorbs+1)/2)
!!  end do
!!end do
!!close(500)
!
!call mma_allocate(dmo1xyz,3,numat,label='dmo1xyz')
!dmo1xyz(1:3,1:numat) = Zero
!
!!---------------------------------------------------
!! partial backtransform one electron density matrix
!
!call mma_allocate(dm1_act,naorbs,naorbs,label='dm1_act')
!dm1_act(1:naorbs,1:naorbs) = Zero
!call density_ci_one(dm1_act)
!
!!write(nf2,*) 'the new transformed dm1'
!
!do i=1,naorbs
!  do j=1,naorbs
!    dm1_act(i,j) = dm1_act(i,j)+Two*p(i,j)
!  end do
!end do
!
!!do i=1,naorbs
!!  do j=1,i
!!    write(nf2,'(2i8,f18.10)') i,j,dm1_act(i,j)
!!  end do
!!end do
!
!do i=1,3
!  do j=1,numat
!    do k=1,naorbs
!      do l=1,k
!        if (k == l) then
!          nkl = ican_a(k)+k
!          dxyz(i,j) = dxyz(i,j)+dsaos(i,j,nkl)*dm1_act(k,l)
!          dmo1xyz(i,j) = dmo1xyz(i,j)+dsaos(i,j,nkl)*dm1_act(k,l)
!
!        else
!          nkl = ican_a(k)+l
!          dxyz(i,j) = dxyz(i,j)+dsaos(i,j,nkl)*dm1_act(k,l)*Two
!
!          dmo1xyz(i,j) = dmo1xyz(i,j)+dsaos(i,j,nkl)*dm1_act(k,l)*Two
!        end if
!      end do
!    end do
!  end do
!end do
!call mma_deallocate(dsaos)
!call mma_deallocate(dm1_act)
!
!call mma_allocate(dgxyz,3,numat,label='dgxyz')
!!dgxyz(1:3,1:numat) = Zero
!!do i=1,numat
!!  do j=1,3
!!    dgxyz(j,i) = dmo1xyz(j,i)*auTokcalmol
!!  end do
!!end do
!call mma_deallocate(dmo1xyz)
!
!!write(nf2,'(//10x,"cartesian coordinate derivatives",//3x,"number  atom ",5x,"x",12x,"y",12x,"z",/)')
!!write(nf2,*) 'the one electron gradient'
!!do i=1,numat
!!  write(nf2,'(6x,i6,3f13.6)') i,(dgxyz(j,i),j=1,3)
!!end do
!
!dgxyz(1:3,1:numat) = Zero
!do i=1,numat
!  do j=1,3
!    dgxyz(j,i) = dxyz(j,i)*auTokcalmol
!  end do
!end do
!
!write(u6,1000)
!
!do i=1,numat
!  write(u6,'(6x,i6,3f13.6)') i,(dgxyz(j,i),j=1,3)
!end do
!call mma_deallocate(dgxyz)
!
!return
!
!1000 format(//10x,'cartesian coordinate derivatives',//3x,'number  atom ',5x,'x',12x,'y',12x,'z',/)
!
!end subroutine grad_one_ao

!subroutine grad_one_mo()
!
!use gugaci_global, only: cf, dm1, dxyz, ican_a, naorbs, norb_all, numat
!use stdalloc, only: mma_allocate, mma_deallocate
!use Constants, only: Zero, Two, auTokcalmol
!use Definitions, only: wp, iwp, u6
!
!implicit none
!integer(kind=iwp) :: i, i0, j, j0, k, l, nkl
!real(kind=wp) :: val
!real(kind=wp), allocatable :: dgxyz(:,:), dmo1xyz(:,:), dsaos(:,:,:)
!
!return
!
!call mma_allocate(dsaos,3,numat,naorbs*(naorbs+1)/2,label='dsaos')
!dsaos(1:3,1:numat,1:naorbs*(naorbs+1)/2) = Zero
!
!!open(500,file='dfock1',form='unformatted')
!!do i=1,numat
!!  do k=1,3
!!    do j=1,naorbs*(naorbs+1)/2
!!      read(500) dsaos(k,i,j)
!!    end do
!!  end do
!!end do
!!close(500)
!
!call mma_allocate(dmo1xyz,3,numat,label='dmo1xyz')
!dmo1xyz(1:3,1:numat) = Zero
!
!do i=1,3
!  do j=1,numat
!
!    do i0=1,norb_all
!      do j0=1,i0
!
!        val = Zero
!        do k=1,naorbs
!          do l=1,k
!            nkl = ican_a(k)+l
!
!            if (k == l) then
!              val = val+dsaos(i,j,nkl)*cf(k,i0)*cf(l,j0)
!            else
!              val = val+dsaos(i,j,nkl)*cf(k,i0)*cf(l,j0)+dsaos(i,j,nkl)*cf(k,j0)*cf(l,i0)
!            end if
!          end do
!        end do
!        if (i0 == j0) then
!          dxyz(i,j) = dxyz(i,j)+dm1(i0,j0)*val
!          dmo1xyz(i,j) = dmo1xyz(i,j)+dm1(i0,j0)*val
!
!        else
!          dxyz(i,j) = dxyz(i,j)+dm1(i0,j0)*val*Two
!          dmo1xyz(i,j) = dmo1xyz(i,j)+dm1(i0,j0)*val*Two
!
!        end if
!        write(u6,'(2i4,2f18.10)') i0,j0,dm1(i0,j0),val
!      end do
!    end do
!  end do
!end do
!call mma_deallocate(dsaos)
!
!call mma_allocate(dgxyz,3,numat,label='dgxyz')
!dgxyz(1:3,1:numat) = Zero
!do i=1,numat
!  do j=1,3
!    dgxyz(j,i) = dmo1xyz(j,i)*auTokcalmol
!  end do
!end do
!call mma_deallocate(dmo1xyz)
!
!write(u6,1000)
!
!do i=1,numat
!  write(u6,'(6x,i6,3f13.6)') i,(dgxyz(j,i),j=1,3)
!end do
!
!dgxyz(1:3,1:numat) = Zero
!do i=1,numat
!  do j=1,3
!    dgxyz(j,i) = dxyz(j,i)*auTokcalmol
!  end do
!end do
!
!write(u6,1000)
!
!do i=1,numat
!  write(u6,'(6x,i6,3f13.6)') i,(dgxyz(j,i),j=1,3)
!end do
!call mma_deallocate(dgxyz)
!
!return
!
!1000 format(//10x,'cartesian coordinate derivatives',//3x,'number  atom ',5x,'x',12x,'y',12x,'z',/)
!
!end subroutine grad_one_mo

!subroutine density_ci_one(dm1_act)
!
!use gugaci_global, only: cf, dm1, naorbs, norb_all, norb_frz
!use Constants, only: Zero
!use Definitions, only: wp, iwp
!
!implicit none
!real(kind=wp), intent(out) :: dm1_act(naorbs,naorbs)
!integer(kind=iwp) :: i, j, norbf, np, nq
!
!norbf = norb_frz+1
!do i=1,naorbs
!  do j=1,i
!    dm1_act(i,j) = Zero
!    do np=norbf,norb_all
!      do nq=norbf,np
!        if (np == nq) then
!          dm1_act(i,j) = dm1_act(i,j)+dm1(np,nq)*cf(i,np)*cf(j,nq)
!        else
!          dm1_act(i,j) = dm1_act(i,j)+dm1(np,nq)*cf(i,np)*cf(j,nq)+dm1(np,nq)*cf(j,np)*cf(i,nq)
!        end if
!      end do
!    end do
!    dm1_act(j,i) = dm1_act(i,j)
!  end do
!end do
!
!end subroutine density_ci_one

!subroutine lagran_fock(x1e,fock)
!
!use gugaci_global, only: ican_a, ican_b, vector1
!use Constants, only: Zero, Two
!use Definitions, only: wp, iwp
!
!implicit none
!real(kind=wp), intent(in) :: x1e(50000)
!real(kind=wp), intent(out) :: fock(n_all,n_all)
!integer(kind=iwp) :: i, i0, j, j0, k, k0, nij, nijkk, nik, nikjk, njk, nkk
!real(kind=wp) :: val
!
!fock(1:n_all,1:n_all) = Zero
!
!do i=1,n_all
!  do j=1,i
!    nij = ican_a(i)+j
!    fock(i,j) = x1e(nij)
!    val = Zero
!    do k=1,n_frz
!      nkk = ican_a(k)+k
!      if (nij >= nkk) then
!        nijkk = ican_b(nij)+nkk
!      else
!        nijkk = ican_b(nkk)+nij
!      end if
!      val = val+vector1(nijkk)*Two
!
!      i0 = max(i,k)
!      k0 = min(i,k)
!      nik = ican_a(i0)+k0
!      j0 = max(j,k)
!      k0 = min(j,k)
!      njk = ican_a(j0)+k0
!      nikjk = ican_b(nik)+njk
!      val = val-vector1(nikjk)
!    end do
!    fock(i,j) = fock(i,j)+val
!    fock(j,i) = fock(i,j)
!  end do
!end do
!
!end subroutine lagran_fock
