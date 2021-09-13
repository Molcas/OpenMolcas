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

subroutine ci_grad()

implicit none
real*8, parameter :: htoklm = 627.50956d0, zero = 0.0d0
#include "drt_h.fh"
#include "intsort_h.fh"
#include "pl_structure_h.fh"
#include "grad_h.fh"
#include "scratch.fh"
#include "lgrn.fh"
#include "iaib.fh"
#include "vect.fh"
#include "grad_xyz.fh"
#include "ncprhf.fh"

!================================================
! the main subroutine for ci gradient calculations.
! ican_a and ican_b save canonical order.
! e(*) save the scf orbital energies.
! xlgrn(*) save the lagrangian matrix.

! ci gradient is still in development, so we do not include it in current
! version of xian-ci

return

!write(6,'(a18,2x,f10.2,2x,a1)') 'end of grad, takes',sc5-sc0,'s'
end subroutine ci_grad

subroutine convert_vector()

implicit none
#include "drt_h.fh"

!=====================================================
! this just uses at debug.

!open(10,file='cigmsvector')
!do i=1,nci_dim
!  read(10,*) iiii,val1,j,val2
!  vector1(j) = val1
!  if ((val2 > 0.0d0) .and. (val1 < 0.0d0)) vector1(j) = -vector1(j)
!  if ((val2 < 0.0d0) .and. (val1 > 0.0d0)) vector1(j) = -vector1(j)
!end do
!close(10)

end subroutine convert_vector

subroutine moread(ii,jj,kk,ll,val)

implicit none
integer :: ii, jj, kk, ll
real*8 :: val
integer :: lri, lrj, lrk, lrl, lrn, nij, nijkl, nkl
#include "drt_h.fh"
#include "grad_h.fh"
#include "iaib.fh"
!====================================================
! transfer the ci mo indices to scf mo indices
! and use library function ishft() to save
! two electron mo indices to dm2_index().

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
nijkl = ican_b(nij)+nkl
vector1(nijkl) = val
!write(nf2,'(4i4,f18.10)') lri,lrj,lrk,lrl,vector1(nijkl)
end subroutine moread

subroutine trans_ijkl_intpos(ii,jj,kk,ll,nxo)

implicit none
integer :: ii, jj, kk, ll, nxo
integer :: lri, lrj, lrk, lrl, lrn, nij, nkl
#include "drt_h.fh"
#include "grad_h.fh"
#include "iaib.fh"
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

subroutine trans_intpos_ijkl(intpos,lijkl)

implicit none
integer :: intpos, lijkl(4)
#include "drt_h.fh"
#include "grad_h.fh"

!=========================================================
! read the saved two electron mo indices.

return
! Avoid unused argument warnings
if (.false.) then
  call Unused_integer(intpos)
  call Unused_integer_array(lijkl)
end if

end subroutine trans_intpos_ijkl

subroutine lagran_act(x1e)

implicit none
#include "drt_h.fh"
#include "grad_h.fh"
#include "vect.fh"
#include "ncprhf.fh"
#include "iaib.fh"
#include "lgrn.fh"
real*8 :: x1e(50000)
integer :: i, i0, j, j0, k, k0, kl, l, l0, m, mik, mjk, nik, nil, niljk, nimkl, nji, njikl, njk, njmkl, nkl, norbf
real*8 :: dum, dumtmp, fock(n_all,n_all)
real*8, parameter :: two = 2.0d0, zero = 0.0d0

!================================================
! lyb
! xlgrn(norb_all,norb_all) is the lagrange matrix.

norbf = n_frz+1
!norbf = 1

xlgrn(1:n_all,1:n_all) = 0.0d0

call lagran_fock(x1e,fock)

! form two electron contributions to the lagrangian with frozen mo

do i=ndbl+1,n_all
  do j=1,n_frz
    dum = fock(j,i)*two
    do k=norbf,n_all
      do l=norbf,n_all
        i0 = max(i,j)
        j0 = min(i,j)
        nji = ican_a(i0)+j0
        k0 = max(k,l)
        l0 = min(k,l)
        nkl = ican_a(k0)+l0
        if (nji >= nkl) then
          njikl = ican_b(nji)+nkl
        else
          njikl = ican_b(nkl)+nji
        end if
        dum = dum+dm1(k0,l0)*vector1(njikl)*two

        i0 = max(i,l)
        l0 = min(i,l)
        nil = ican_a(i0)+l0
        j0 = max(j,k)
        k0 = min(j,k)
        njk = ican_a(j0)+k0
        if (nil >= njk) then
          niljk = ican_b(nil)+njk
        else
          niljk = ican_b(njk)+nil
        end if
        dum = dum-dm1(k0,l0)*vector1(niljk)
      end do
    end do
    xlgrn(i,j) = xlgrn(i,j)+dum
  end do
end do

do i=1,n_frz
  do j=ndbl+1,n_all
    dum = fock(j,i)*two
    do k=norbf,n_all
      do l=norbf,n_all
        i0 = max(i,j)
        j0 = min(i,j)
        nji = ican_a(i0)+j0
        k0 = max(k,l)
        l0 = min(k,l)
        nkl = ican_a(k0)+l0
        if (nji >= nkl) then
          njikl = ican_b(nji)+nkl
        else
          njikl = ican_b(nkl)+nji
        end if
        dum = dum+dm1(k0,l0)*vector1(njikl)*two

        i0 = max(i,l)
        l0 = min(i,l)
        nil = ican_a(i0)+l0
        j0 = max(j,k)
        k0 = min(j,k)
        njk = ican_a(j0)+k0
        if (nil >= njk) then
          niljk = ican_b(nil)+njk
        else
          niljk = ican_b(njk)+nil
        end if
        dum = dum-dm1(k0,l0)*vector1(niljk)
      end do
    end do
    xlgrn(i,j) = xlgrn(i,j)+dum
  end do
end do

! form two electron contributions to the lagrangian with active mo

do i=norbf,ndbl
  do j=ndbl+1,n_all
    dum = zero
    do m=norbf,n_all
      i0 = ican_a(max(i,m))+min(i,m)
      j0 = ican_a(max(j,m))+min(j,m)
      kl = 0
      do k=norbf,n_all
        dumtmp = zero
        do l=norbf,k-1
          kl = kl+1
          if (kl > i0) then
            nimkl = ican_b(kl)+i0
          else
            nimkl = ican_b(i0)+kl
          end if

          if (kl > j0) then
            njmkl = ican_b(kl)+j0
          else
            njmkl = ican_b(j0)+kl
          end if
          dumtmp = dumtmp+vector1(nimkl)*vector2(njmkl)
          !if ((i == 1) .and. (j == 5)) write(nf2,'(6i4,i8,f18.10)') i,j,m,k,l,kl,njmkl,vector2(njmkl)
        end do

        dum = dum+dumtmp*two

        kl = kl+1
        if (kl > i0) then
          nimkl = ican_b(kl)+i0
        else
          nimkl = ican_b(i0)+kl
        end if

        if (kl > j0) then
          njmkl = ican_b(kl)+j0
        else
          njmkl = ican_b(j0)+kl
        end if

        dum = dum+vector1(nimkl)*vector2(njmkl)
        !if (nimkl /= njmkl) write(nf2,'(2i8,2f18.10)') nimkl,njmkl,vector1(nimkl),vector2(njmkl)
        !if ((i == 1) .and. (j == 5)) write(nf2,'(6i4,i8,f18.10)') i,j,m,k,l,kl,njmkl,vector2(njmkl
      end do
    end do

    xlgrn(i,j) = xlgrn(i,j)+dum*two
    !write(2,'(a7,2i4,f18.10)') 'xlgrn_2',i,j,xlgrn(i,j)
  end do
end do
do i=norbf,ndbl
  do j=ndbl+1,n_all
    dum = zero
    do k=norbf,norb_all
      mik = max(i,k)
      nik = min(i,k)
      mjk = max(j,k)
      njk = min(j,k)
      dum = dum+fock(mik,nik)*dm1(mjk,njk)

      !write(nf2,'(2f18.10)') x1e(mnik),dm1(mjk,njk)
    end do
    xlgrn(i,j) = xlgrn(i,j)+dum
    !write(2,'(a9,2i4,f18.10)') 'xlgrn_all',i,j,xlgrn(i,j)
  end do
end do

do i=ndbl+1,n_all
  do j=norbf,ndbl
    dum = zero
    do m=norbf,n_all
      i0 = ican_a(max(i,m))+min(i,m)
      j0 = ican_a(max(j,m))+min(j,m)
      kl = 0
      do k=norbf,n_all
        dumtmp = zero
        do l=norbf,k-1
          kl = kl+1
          if (kl > i0) then
            nimkl = ican_b(kl)+i0
          else
            nimkl = ican_b(i0)+kl
          end if

          if (kl > j0) then
            njmkl = ican_b(kl)+j0
          else
            njmkl = ican_b(j0)+kl
          end if
          dumtmp = dumtmp+vector1(nimkl)*vector2(njmkl)
          !if ((i == 1) .and. (j == 5)) write(nf2,'(6i4,i8,f18.10)') i,j,m,k,l,kl,njmkl,vector2(njmkl)
        end do

        dum = dum+dumtmp*two

        kl = kl+1
        if (kl > i0) then
          nimkl = ican_b(kl)+i0
        else
          nimkl = ican_b(i0)+kl
        end if

        if (kl > j0) then
          njmkl = ican_b(kl)+j0
        else
          njmkl = ican_b(j0)+kl
        end if

        dum = dum+vector1(nimkl)*vector2(njmkl)
        !if (nimkl /= njmkl) write(nf2,'(2i8,2f18.10)') nimkl,njmkl,vector1(nimkl),vector2(njmkl)
        !if ((i == 1) .and. (j == 5)) write(nf2,'(6i4,i8,f18.10)') i,j,m,k,l,kl,njmkl,vector2(njmkl
      end do
    end do

    xlgrn(i,j) = xlgrn(i,j)+dum*two
    !write(2,'(a7,2i4,f18.10)') 'xlgrn_2',i,j,xlgrn(i,j)
  end do
end do

do i=ndbl+1,n_all
  do j=norbf,ndbl
    dum = zero
    do k=norbf,norb_all
      mik = max(i,k)
      nik = min(i,k)
      mjk = max(j,k)
      njk = min(j,k)
      dum = dum+fock(mik,nik)*dm1(mjk,njk)

      !write(nf2,'(2f18.10)') x1e(mnik),dm1(mjk,njk)
    end do
    xlgrn(i,j) = xlgrn(i,j)+dum
    !write(2,'(a9,2i4,f18.10)') 'xlgrn_all',i,j,xlgrn(i,j)
  end do
end do

end subroutine lagran_act

subroutine lagran_all(x1e)

implicit none
real*8 :: x1e(50000)
integer :: i, i0, j, j0, k, kl, l, m, mik, mjk, mnik, nik, nimkl, njk, njmkl
real*8 :: dum, dumtmp
real*8, parameter :: two = 2.0d0, zero = 0.0d0
#include "drt_h.fh"
#include "grad_h.fh"
#include "iaib.fh"
#include "lgrn.fh"
#include "ncprhf.fh"

!================================================
! lyb
! xlgrn(norb_all,norb_all) is the lagrange matrix.

!norbf = n_frz+1
!norbf = 1

xlgrn(1:norb_all,1:norb_all) = 0.0d0

! form two electron contributions to the lagrangian
do i=1,norb_all
  do j=1,norb_all
    dum = zero
    do m=1,norb_all
      i0 = ican_a(max(i,m))+min(i,m)
      j0 = ican_a(max(j,m))+min(j,m)
      kl = 0
      do k=1,norb_all
        dumtmp = zero
        do l=1,k-1
          kl = kl+1
          if (kl > i0) then
            nimkl = ican_b(kl)+i0
          else
            nimkl = ican_b(i0)+kl
          end if

          if (kl > j0) then
            njmkl = ican_b(kl)+j0
          else
            njmkl = ican_b(j0)+kl
          end if
          dumtmp = dumtmp+vector1(nimkl)*vector2(njmkl)
          !if ((i == 1) .and. (j == 5)) write(nf2,'(6i4,i8,f18.10)') i,j,m,k,l,kl,njmkl,vector2(njmkl)
        end do

        dum = dum+dumtmp*two

        kl = kl+1
        if (kl > i0) then
          nimkl = ican_b(kl)+i0
        else
          nimkl = ican_b(i0)+kl
        end if

        if (kl > j0) then
          njmkl = ican_b(kl)+j0
        else
          njmkl = ican_b(j0)+kl
        end if

        dum = dum+vector1(nimkl)*vector2(njmkl)
        !if (nimkl /= njmkl) write(nf2,'(2i8,2f18.10)') nimkl,njmkl,vector1(nimkl),vector2(njmkl)
        !if ((i == 1) .and. (j == 5)) write(nf2,'(6i4,i8,f18.10)') i,j,m,k,l,kl,njmkl,vector2(njmkl)
      end do
    end do

    xlgrn(i,j) = xlgrn(i,j)+dum*two
    !write(2,'(a7,2i4,f18.10)') 'xlgrn_2',i,j,xlgrn(i,j)
  end do
end do
do i=1,norb_all
  do j=1,norb_all
    dum = zero
    do k=1,norb_all
      mik = max(i,k)
      nik = min(i,k)
      mjk = max(j,k)
      njk = min(j,k)
      mnik = ican_a(mik)+nik
      dum = dum+x1e(mnik)*dm1(mjk,njk)

      !write(nf2,'(2f18.10)') x1e(mnik),dm1(mjk,njk)
    end do
    xlgrn(i,j) = xlgrn(i,j)+dum
    !write(2,'(a9,2i4,f18.10)') 'xlgrn_all',i,j,xlgrn(i,j)
  end do
end do

end subroutine lagran_all

subroutine writedm2(nx)

implicit none
integer :: nx
!character(len=256) :: filename
#include "drt_h.fh"
#include "grad_h.fh"
#include "scratch.fh"

!filename = tmpdir(1:len_str)//'/density'
!len = len_str+8
! Need debug
!open(nf22,file=filename(1:len),form='unformatted')

!open(nf22,file='density',form='unformatted')
!write(nf22) (vector2(i),i=1,nx)
!close(nf22)

! Avoid unused argument warnings
if (.false.) call Unused_integer(nx)

end subroutine writedm2

subroutine readdm2(nx)

implicit none
integer :: nx
!character(len=256) :: filename
#include "drt_h.fh"
#include "grad_h.fh"
#include "scratch.fh"

vector2(1:nx) = 0.0d+00

!filename = tmpdir(1:len_str)//'/density'
!len = len_str+8

!open(nf22,file=filename(1:len),form='unformatted')

!open(nf22,file='density',form='unformatted')
!read(nf22) (vector2(i),i=1,nx)
!close(nf22)

end subroutine readdm2

#ifdef _COMPILE_

subroutine backtransmo()

implicit none
#include "drt_h.fh"
#include "grad_h.fh"
#include "scratch.fh"
#include "vect.fh"
#include "iaib.fh"
#include "density.fh"
integer :: ican_ab(norb_all)
real*8 :: c(70000), dm1_act(naorbs,naorbs)
logical :: resina
!character(len=256) :: filename
real*8, parameter :: half = 0.5d0, half2 = 0.25d0, htokl = 627.50956d0, zero = 0.0d0, two = 2.0d0

!***********************************************************************
!
!   program transmo transfers the one and two electronic ao integrals to
!   integrals.
!          on input
! cf      = all m.o.'s coefficient
! naorbs  = number of atomic orbitals
! maxao   = the number of ao integrals
! h       = the one electronic ao integrals
! vector1 = the two electronic ao density matrix
! norb    = the start mo wanted to be transfered
! norbs   = the end mo wanted to be transfered
! c       = the m.o.'s coefficient with the order that is different from
!           this is need by the transformation subroutine c4itd
!
!          on output
! tmoint1 = the one electronic mo integrals
! vector2 = the two electronic mo density matrix
!
! c4itd   = the subroutine to transfer aos to mos
!           copy right by carlos f. bunge. annik vivier bunge. gerardo c
!           and jean-pierre daudey, 1987.
!           reference: c. f. bunge, a. v. bunge, g. cisneros and j.-p. d
!                      comput. chem. vol 12, page 91, year 1988.
! iteifd  = the subroutine to scale the aos by a factor 1/2
!***********************************************************************

time0 = c_time()

!write(nf2,*) 'start of backtransform'

norbf = norb_frz

norbs1 = naorbs*(naorbs+1)*0.5
norbe = norb_all-norbf
norbs2 = norbe*(norbe+1)*0.5
norbs3 = (norbe-norbf+1)*naorbs
nx = (naorbs*naorbs+naorbs)*(naorbs*naorbs+naorbs+2)*0.125
mx = (norbe*norbe+norbe)*(norbe*norbe+norbe+2)*0.125

vector1(1:nx) = zero

l1 = norbe
do i=1,l1
  i0 = i-1
  ican_ab(i) = i0*l1-(i0*i0-i0)/2
end do

ij = 0
do i=norbf+1,norb_all
  do j=1,naorbs
    ij = ij+1
    c(ij) = cf(j,i)
    !write(nf2,'(2i4,2f18.10)') i,j,c(ij),cf(j,i)
  end do
end do

do i0=norbf+1,norb_all
  do j0=norbf+1,i0
    do k0=norbf+1,i0
      l0max = k0
      if (i0 == k0) l0max = j0
      do l0=norbf+1,l0max
        if (l0 <= norbf) cycle
        nij0 = ican_a(i0)+j0
        nkl0 = ican_a(k0)+l0
        nijkl0 = ican_b(nij0)+nkl0
        val = vector2(nijkl0)
        i = i0-norbf
        j = j0-norbf
        k = k0-norbf
        l = l0-norbf

        nij = ican_ab(j)+i-j+1
        nkl = ican_ab(l)+k-l+1
        if (nij >= nkl) then
          nijkl = ican_b(nij)+nkl
        else
          nijkl = ican_b(nkl)+nij
        end if
        vector1(nijkl) = val
        !write(nf2,'(8i4,2i8,f18.10)') i0,j0,k0,l0,i,j,k,l,nijkl0,nijkl,val

      end do
    end do
  end do
end do

!***********************************************************************
! the aoints should be saved as the order the following
! provided.
!
! nijkl = 0
! do i=1,naorbs
!   do j=i,naorbs
!     do k=1,i
!       if (k == i) then
!         inl = j
!       else
!         inl = naorbs
!       end if
!       do l=k,inl
!         nijkl = nijkl+1
!         read(naoint) vector1(nijkl)
!
!         write(6,'(4i3,i8,f18.10)') i,j,k,l,nijkl,vector1(nijkl)
!         write(nf2,'(5i8)') i,j,k,l,nijkl
!
!       end do
!     end do
!   end do
! end do
!***********************************************************************

vector2(1:nx) = zero
resina = .false.
nsym = 1
ncase = 1

call iteifd(ncase,nsym,norbe,norbe,norbe,norbe,vector1)

call c4itd(norbe,norbe,norbe,norbe,naorbs,naorbs,naorbs,naorbs,nsym,ncase,c,c,c,c,vector1,resina,vector2)

vector2(1:nx) = zero

!-----------------------------------------------------------------------
dm1_act(1:naorbs,1:naorbs) = zero

call density_ci_one(dm1_act)

!do i=1,naorbs
!  do j=1,naorbs
!    write(nf2,'(2i4,f18.10)') i,j,dm1_act(i,j)
!  end do
!end do

nijkl = 0
do i=1,naorbs
  do j=1,i
    do k=1,i
      if (k == i) then
        inl = j
      else
        inl = k
      end if
      do l=1,inl
        nijkl = nijkl+1

        valtmp = two*p(i,j)*p(k,l)-half*p(i,l)*p(j,k)-half*p(i,k)*p(j,l)+p(i,j)*dm1_act(k,l)+p(k,l)*dm1_act(i,j)- &
                 half2*(p(j,k)*dm1_act(i,l)+p(j,l)*dm1_act(i,k)+p(i,k)*dm1_act(j,l)+p(i,l)*dm1_act(j,k))

        vector1(nijkl) = vector1(nijkl)+valtmp
        !write(nf2,'(4i3,i8,f18.10)') i,j,k,l,nijkl,vector1(nijkl)
        !write(nf2,'(5i8)') i,j,k,l,nijkl

      end do
    end do
  end do
end do

!write(nf2,*) 'the new dm2'

!do i=1,nx
!  write(nf2,'(i8,f18.10)')i, vector1(i)
!end do

!filename = tmpdir(1:len_str)//'/backdm2'
!len = len_str+8

!open(20,file=filename(1:len),form='unformatted')

!open(20,file='backdm2',form='unformatted')
!write(20) (vector1(i),i=1,nx)
!close(20)

time1 = c_time()-time0
!write(nf2,'(4x,"trans run time =",f8.3,2x,"seconds")') time1
!write(nf2,*) 'end of backtransform'

100 format(2i4,f18.10)
200 format(4i4,4x,f18.10)
300 format(4x,a6,2x,f18.10,5i4,f20.15)

end subroutine backtransmo

#endif

subroutine backtrans_test()

implicit none
integer :: i, i0, j, j0, k, k0, l, l0, nij, nijkl, nkl
real*8 :: sum_, val
#include "drt_h.fh"
#include "vect.fh"
#include "iaib.fh"

!-------------------------------------------------------------
! this subroutine for test the backtrans result by just one density matr
! ao, such as (ij|kl)

sum_ = 0.0d+00
i0 = 6
j0 = 6
k0 = 6
l0 = 6

do i=norb_frz+1,norb_all
  do j=norb_frz+1,norb_all
    if (i >= j) then
      nij = ican_a(i)+j
    else
      nij = ican_a(j)+i
    end if

    do k=norb_frz+1,norb_all
      do l=norb_frz+1,norb_all
        if (k >= l) then
          nkl = ican_a(k)+l
        else
          nkl = ican_a(l)+k
        end if
        if (nij >= nkl) then
          nijkl = ican_b(nij)+nkl
        else
          nijkl = ican_b(nkl)+nij
        end if
        val = vector2(nijkl)
        sum_ = sum_+val*cf(i0,i)*cf(j0,j)*cf(k0,k)*cf(l0,l)
      end do
    end do
  end do
end do

!write(nf2,*) 'num 1 sum=',sum_

end subroutine backtrans_test

! FIXME: ndao is undefined
!subroutine grad_two()
!
!implicit none
!#include "drt_h.fh"
!#include "grad_xyz.fh"
!#include "iaib.fh"
!#include "scratch.fh"
!integer :: index_atom(3,numat*(numat+1)/2), ndi0(ndao), ndj0(ndao), ndk0(ndao), ndl0(ndao)
!real*8 :: daoint1(ndao), daoxyz(3,numat), dgxyz(3,numat)
!!character(len=256) filename
!real*8, parameter :: four=4.0d0, htoklm = 627.50956d0, one=1.0d0, two=2.0d0, zero=0.0d0
!
!npat = numat*(numat+1)/2
!index_atom(1:3,1:npat) = 0
!
!ndi0(:) = 0
!ndj0(:) = 0
!ndk0(:) = 0
!ndl0(:) = 0
!daoint1(:) = zero
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
!do i=1,numat
!  do j=1,i-1
!    nnij = ican_a(i)+j
!    !write(nf2,*) i,ican_a(i)
!    do k=1,3
!      ind = index_atom(k,nnij)
!      val = zero
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
!        aa = one
!        bb = one
!        if (i0 /= j0) aa = two*aa
!        if (k0 /= l0) bb = two*bb
!        !===================================================
!        ! this place should multiple two, because the gradient
!        ! integral always exist the relation that :
!        ! index nij=i*(i-1)/2+j
!        ! index nkl=k*(k-1)/2+l
!        ! but here always the nij >= nkl.
!
!        val = val+val1*val2*aa*bb*two
!      end do
!      daoxyz(k,i) = daoxyz(k,i)+val
!      daoxyz(k,j) = daoxyz(k,j)-val
!
!      dxyz(k,i) = dxyz(k,i)+val
!      dxyz(k,j) = dxyz(k,j)-val
!    end do
!  end do
!end do
!
!!dgxyz(1:3,1:numat) = zero
!!do i=1,numat
!!  do j=1,3
!!    dgxyz(j,i) = daoxyz(j,i)*htoklm
!!  end do
!!end do
!
!
!!write(nf2,'(//10x,'cartesian coordinate derivatives',//3x,'number  atom ',5x,'x',12x,'y',12x,'z',/)')
!
!!do i=1,numat
!!  write(nf2,'(6x,i6,3f13.6)') i,(dgxyz(j,i),j=1,3)
!!end do
!
!dgxyz(1:3,1:numat) = zero
!do i=1,numat
!  do j=1,3
!    dgxyz(j,i) = dxyz(j,i)*htoklm
!  end do
!end do
!
!
!!write(nf2,'(//10x,'cartesian coordinate derivatives',//3x,'number  atom ',5x,'x',12x,'y',12x,'z',/)')
!
!do i=1,numat
!  write(6,'(6x,i6,3f13.6)') i,(dgxyz(j,i),j=1,3)
!end do
!
!end subroutine backtrans_test

subroutine grad_one_ao()

implicit none
#include "drt_h.fh"
#include "grad_h.fh"
#include "scratch.fh"
#include "lgrn.fh"
#include "iaib.fh"
#include "vect.fh"
#include "grad_xyz.fh"
#include "density.fh"
integer :: i, j, k, l, nkl
real*8 :: dgxyz(3,numat), dm1_act(naorbs,naorbs), dmo1xyz(3,numat), dsaos(3,numat,naorbs*(naorbs+1)/2)
!character*256 filename
real*8, parameter :: four = 4.0d0, htoklm = 627.50956d0, one = 1.0d0, two = 2.0d0, zero = 0.0d0

dsaos(1:3,1:numat,1:naorbs*(naorbs+1)/2) = zero

!filename = tmpdir(1:len_str)//'/dfock1'
!len = len_str+7

!open(500,file=filename(1:len),form='unformatted')

!open(500,file='dfock1',form='unformatted')
!do i=1,numat
!  do k=1,3
!    read(500) (dsaos(k,i,j),j=1,naorbs*(naorbs+1)/2)
!  end do
!end do
!close(500)

dmo1xyz(1:3,1:numat) = zero

!---------------------------------------------------
! partial backtransform one electron density matrix

dm1_act(1:naorbs,1:naorbs) = zero
call density_ci_one(dm1_act)

!write(nf2,*) 'the new transformed dm1'

do i=1,naorbs
  do j=1,naorbs
    dm1_act(i,j) = dm1_act(i,j)+two*p(i,j)
  end do
end do

!do i=1,naorbs
!  do j=1,i
!    write(nf2,'(2i8,f18.10)') i,j,dm1_act(i,j)
!  end do
!end do

do i=1,3
  do j=1,numat
    do k=1,naorbs
      do l=1,k
        if (k == l) then
          nkl = ican_a(k)+k
          dxyz(i,j) = dxyz(i,j)+dsaos(i,j,nkl)*dm1_act(k,l)
          dmo1xyz(i,j) = dmo1xyz(i,j)+dsaos(i,j,nkl)*dm1_act(k,l)

        else
          nkl = ican_a(k)+l
          dxyz(i,j) = dxyz(i,j)+dsaos(i,j,nkl)*dm1_act(k,l)*two

          dmo1xyz(i,j) = dmo1xyz(i,j)+dsaos(i,j,nkl)*dm1_act(k,l)*two
        end if
      end do
    end do
  end do
end do

!dgxyz(1:3,1:numat)  =zero
!do i=1,numat
!  do j=1,3
!    dgxyz(j,i) = dmo1xyz(j,i)*htoklm
!  end do
!end do

!write(nf2,'(//10x,"cartesian coordinate derivatives",//3x,"number  atom ",5x,"x",12x,"y",12x,"z",/)')
!write(nf2,*) 'the one electron gradient'
!do i=1,numat
!  write(nf2,'(6x,i6,3f13.6)') i,(dgxyz(j,i),j=1,3)
!end do

dgxyz(1:3,1:numat) = zero
do i=1,numat
  do j=1,3
    dgxyz(j,i) = dxyz(j,i)*htoklm
  end do
end do

write(6,1000)

do i=1,numat
  write(6,'(6x,i6,3f13.6)') i,(dgxyz(j,i),j=1,3)
end do

1000 format(//10x,'cartesian coordinate derivatives',//3x,'number  atom ',5x,'x',12x,'y',12x,'z',/)

end subroutine grad_one_ao

subroutine grad_one_mo()

implicit none
#include "drt_h.fh"
#include "grad_h.fh"
#include "lgrn.fh"
#include "iaib.fh"
#include "vect.fh"
#include "grad_xyz.fh"
integer :: i, i0, j, j0, k, l, nkl
real*8 :: dgxyz(3,numat), dmo1xyz(3,numat), dsaos(3,numat,naorbs*(naorbs+1)/2), val
real*8, parameter :: four = 4.0d0, htoklm = 627.50956d0, one = 1.0d0, two = 2.0d0, zero = 0.0d0

return

dsaos(1:3,1:numat,1:naorbs*(naorbs+1)/2) = zero

!open(500,file='dfock1',form='unformatted')
!do i=1,numat
!  do k=1,3
!    do j=1,naorbs*(naorbs+1)/2
!      read(500) dsaos(k,i,j)
!    end do
!  end do
!end do
!close(500)

dmo1xyz(1:3,1:numat) = zero

do i=1,3
  do j=1,numat

    do i0=1,norb_all
      do j0=1,i0

        val = zero
        do k=1,naorbs
          do l=1,k
            nkl = ican_a(k)+l

            if (k == l) then
              val = val+dsaos(i,j,nkl)*cf(k,i0)*cf(l,j0)
            else
              val = val+dsaos(i,j,nkl)*cf(k,i0)*cf(l,j0)+dsaos(i,j,nkl)*cf(k,j0)*cf(l,i0)
            end if
          end do
        end do
        if (i0 == j0) then
          dxyz(i,j) = dxyz(i,j)+dm1(i0,j0)*val
          dmo1xyz(i,j) = dmo1xyz(i,j)+dm1(i0,j0)*val

        else
          dxyz(i,j) = dxyz(i,j)+dm1(i0,j0)*val*two
          dmo1xyz(i,j) = dmo1xyz(i,j)+dm1(i0,j0)*val*two

        end if
        write(6,'(2i4,2f18.10)') i0,j0,dm1(i0,j0),val
      end do
    end do
  end do
end do

dgxyz(1:3,1:numat) = zero
do i=1,numat
  do j=1,3
    dgxyz(j,i) = dmo1xyz(j,i)*htoklm
  end do
end do

write(6,1000)

do i=1,numat
  write(6,'(6x,i6,3f13.6)') i,(dgxyz(j,i),j=1,3)
end do

dgxyz(1:3,1:numat) = zero
do i=1,numat
  do j=1,3
    dgxyz(j,i) = dxyz(j,i)*htoklm
  end do
end do

write(6,1000)

do i=1,numat
  write(6,'(6x,i6,3f13.6)') i,(dgxyz(j,i),j=1,3)
end do

1000 format(//10x,'cartesian coordinate derivatives',//3x,'number  atom ',5x,'x',12x,'y',12x,'z',/)

end subroutine grad_one_mo

subroutine density_ci_one(dm1_act)

implicit none
#include "drt_h.fh"
#include "grad_h.fh"
#include "iaib.fh"
#include "vect.fh"
integer :: i, j, norbf, np, nq
real*8 :: dm1_act(naorbs,naorbs)
real*8, parameter :: four = 4.0d0, one = 1.0d0, two = 2.0d0, zero = 0.0d0

norbf = norb_frz+1
do i=1,naorbs
  do j=1,i
    dm1_act(i,j) = zero
    do np=norbf,norb_all
      do nq=norbf,np
        if (np == nq) then
          dm1_act(i,j) = dm1_act(i,j)+dm1(np,nq)*cf(i,np)*cf(j,nq)
        else
          dm1_act(i,j) = dm1_act(i,j)+dm1(np,nq)*cf(i,np)*cf(j,nq)+dm1(np,nq)*cf(j,np)*cf(i,nq)
        end if
      end do
    end do
    dm1_act(j,i) = dm1_act(i,j)
  end do
end do

end subroutine density_ci_one

subroutine lagran_fock(x1e,fock)

implicit none
#include "drt_h.fh"
#include "grad_h.fh"
#include "vect.fh"
#include "iaib.fh"
#include "ncprhf.fh"
real*8 :: x1e(50000)
integer :: i, i0, j, j0, k, k0, nij, nijkk, nik, nikjk, njk, nkk
real*8 :: fock(n_all,n_all), val
real*8, parameter :: two = 2.0d0, zero = 0.0d0

fock(1:n_all,1:n_all) = zero

do i=1,n_all
  do j=1,i
    nij = ican_a(i)+j
    fock(i,j) = x1e(nij)
    val = zero
    do k=1,n_frz
      nkk = ican_a(k)+k
      if (nij >= nkk) then
        nijkk = ican_b(nij)+nkk
      else
        nijkk = ican_b(nkk)+nij
      end if
      val = val+vector1(nijkk)*two

      i0 = max(i,k)
      k0 = min(i,k)
      nik = ican_a(i0)+k0
      j0 = max(j,k)
      k0 = min(j,k)
      njk = ican_a(j0)+k0
      nikjk = ican_b(nik)+njk
      val = val-vector1(nikjk)
    end do
    fock(i,j) = fock(i,j)+val
    fock(j,i) = fock(i,j)
  end do
end do

end subroutine lagran_fock
