!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!                                                                      *
! Copyright (C) 2014, Steven Vancoillie                                *
!***********************************************************************

module faroald
! written by Steven Vancoillie, summer 2014
!
! The faroald module handles sigma updates as the result of acting
! with the hamiltonian operator on a CI vector: s = H c, where s and c
! are CI expansions in determinant basis.
!
! The implementation follows the minimum operation count algorithm
! published by Olsen & Co in J. Chem. Phys. 89, 2185 (1988).

use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, Half
use Definitions, only: wp, iwp
#ifdef _PROF_
use, intrinsic :: iso_fortran_env, only: int64
use Definitions, only: u6
#endif

implicit none
private

! wavefunction info
integer(kind=iwp) :: my_nel, my_norb, mult, nela, nelb, nhoa, nhob, ndeta, ndetb, my_ndet

! excitation tables
type ex1_struct
  integer(kind=iwp) :: p, q, sgn, rank
end type

type(ex1_struct), allocatable :: ex1_a(:,:), ex1_b(:,:)
integer(kind=iwp) :: max_ex1a, max_ex1b, max_ex2a, max_ex2b, max_LRs

#ifdef _PROF_
integer(kind=int64) :: nflop
#endif

public :: ex1_a, ex1_b, ex1_init, max_LRs, max_ex1a, max_ex1b, max_ex2a, max_ex2b, mult, my_ndet, my_nel, my_norb, ndeta, ndetb, &
          nela, nelb, nhoa, nhob, sigma_update

! Extensions to mma interfaces

interface cptr2loff
  module procedure ex1_cptr2loff
end interface
interface mma_allocate
  module procedure ex1_mma_allo_2D, ex1_mma_allo_2D_lim
end interface
interface mma_deallocate
  module procedure ex1_mma_free_2D
end interface

public :: mma_allocate, mma_deallocate

contains

subroutine sigma_update(h,g,sgm,psi)
! The sigma update routine performs the following operation:
! |sgm> = H |psi>, with the hamiltonian defined by
! H = sum_tu h(t,u) E_tu + sum_tuvx g(t,u,v,x) E_tuvx.

  ! integrals
  real(kind=wp), intent(in) :: h(my_norb,my_norb), g(my_norb,my_norb,my_norb,my_norb)
  ! wavefunctions
  real(kind=wp), intent(out) :: sgm(:,:)
  real(kind=wp), intent(in) :: psi(:,:)
  integer(kind=iwp) :: t, u, v, & ! orbital indices
                       iasta, iaend, ibsta, ibend ! determinant index ranges
  real(kind=wp), allocatable :: k(:,:), psiT(:,:), sgmT(:,:)
# ifdef _PROF_
  ! profiling
  real(kind=wp) :: t1_cpu, t2_cpu, tot_cpu, t1_wall, t2_wall, tot_wall, walltime, flops
# endif

  ! Distributes a dimension over processes
  ! and returns the range of this process.
  call par_range(ndeta,iasta,iaend)
  call par_range(ndetb,ibsta,ibend)

# ifdef _PROF_
  ! initialize flop count
  nflop = 0

  call timing(t1_cpu,tot_cpu,t1_wall,tot_wall)
# endif

  ! Set the sigma vector to 0
  sgm = Zero

  ! First, construct a new effective one-electron integral matrix:
  ! k_tu = h_tu - 1/2 sum_v g_tvvu, to be used with sigma1/sigma2.
  call mma_allocate(k,my_norb,my_norb,label='k')
  do u=1,my_norb
    do t=1,my_norb
      ! g_tvvu = g_vtvu
      k(t,u) = Zero
      do v=1,my_norb
        k(t,u) = k(t,u)+g(v,t,v,u)
      end do
      k(t,u) = h(t,u)-Half*k(t,u)
    end do
  end do

  ! Second, for sigma2 and sigma3, we are better off with the transpose
  ! of sgm and/or psi, so allocate and assign them here.
  call mma_allocate(psiT,ndetb,ndeta,label='psiT')
  call dtrans(ndeta,ndetb,psi,ndeta,psiT,ndetb)

  ! Now the actual contributions to sigma are computed. For a singlet
  ! (mult = 1), sigma2 is not computed and sigma3 will only do half the
  ! work.  But at the end, the transpose of sigma has to be added to
  ! sigma. This should be more efficient than computing everything, but
  ! note that the transpose operation is also included in the timings
  ! used to compute the flop efficiency.

  call sigma1(k,g,sgm,psi,ibsta,ibend)

  if (mult /= 1) then
    ! we need efficient access to sgm by using the transpose
    call mma_allocate(sgmT,ndetb,ndeta,label='sgmT')
    call dtrans(ndeta,ndetb,sgm,ndeta,sgmT,ndetb)

    call sigma2(k,g,sgmT,psiT,iasta,iaend)

    call dtrans(ndetb,ndeta,sgmT,ndetb,sgm,ndeta)
    call mma_deallocate(sgmT)
  end if

  call sigma3(g,sgm,psiT,ibsta,ibend)

  ! sum over all processes
  call gadsum(sgm,ndeta*ndetb)

  if (mult == 1) then
    ! for Ms = 0 (only used for singlet), sgm := sgm + sgm^T
    call transadd(ndeta,sgm,ndeta)
  end if

  call mma_deallocate(psiT)
  call mma_deallocate(k)

# ifdef _PROF_
  call timing(t2_cpu,tot_cpu,t2_wall,tot_wall)

  walltime = t2_wall-t1_wall

  if (walltime /= Zero) then
    flops = nflop/walltime
    write(u6,'(1x,a,2(f10.3,a))') 'sigma update: ',walltime,' s, ',flops*1.0d-9,' Gflops.'
  end if
# endif

end subroutine sigma_update

subroutine sigma1(k,g,sgm,psi,ibsta,ibend)
! sigma1 = sum_jb sum_tu <jb|E_tu|ib> (h_kl - 1/2 sum_v <tv|vx>) C(ia,jb)
!    + 1/2 sum_jb sum_tuvx <jb|E_tu E_vx|ib> g_tuvx C(ia,jb)

  ! integrals
  real(kind=wp), intent(in) :: k(my_norb,my_norb), g(my_norb,my_norb,my_norb,my_norb)
  ! wavefunctions
  real(kind=wp), intent(inout) :: sgm(:,:)
  real(kind=wp), intent(in) :: psi(:,:)
  integer(kind=iwp), intent(in) :: ibsta, ibend
  ! local variables
  integer(kind=iwp) :: ib, jb, kb, t, u, v, x, tu, vx, sgn_tu, sgn_vx
  real(kind=wp), allocatable :: f(:)

  call mma_allocate(f,ndetb,label='f')

  do ib=ibsta,ibend
    ! f array construction
    f = Zero
    do tu=1,max_ex1b
      t = ex1_b(tu,ib)%p
      u = ex1_b(tu,ib)%q
      sgn_tu = ex1_b(tu,ib)%sgn
      kb = ex1_b(tu,ib)%rank
      f(kb) = f(kb)+sgn_tu*k(t,u)
      do vx=1,max_ex1b
        v = ex1_b(vx,kb)%p
        x = ex1_b(vx,kb)%q
        sgn_vx = ex1_b(vx,kb)%sgn
        jb = ex1_b(vx,kb)%rank
        f(jb) = f(jb)+Half*sgn_tu*sgn_vx*g(v,x,t,u)
      end do
    end do
    ! sigma addition
    kb = 0
    do jb=1,ndetb
      if (f(jb) /= Zero) then
        kb = kb+1
#       ifdef _PROF_
        nflop = nflop+2*ndeta
#       endif
        call daxpy_(ndeta,f(jb),psi(:,jb),1,sgm(:,ib),1)
      end if
    end do
    if (kb > max_ex2b) stop 'exceeded max double excitations'
  end do

  call mma_deallocate(f)

end subroutine sigma1

subroutine sigma2(k,g,sgm,psi,iasta,iaend)
! sigma2 = sum_ja sum_tu <ja|E_tu|ia> (h_kl - 1/2 sum_v <tv|vx>) C(ja,ib)
!    + 1/2 sum_ja sum_tuvx <ja|E_tu E_vx|ia> g_tuvx C(ja,ib)

  ! integrals
  real(kind=wp), intent(in) :: k(my_norb,my_norb), g(my_norb,my_norb,my_norb,my_norb)
  ! wavefunctions
  real(kind=wp), intent(inout) :: sgm(:,:)
  real(kind=wp), intent(in) :: psi(:,:)
  integer(kind=iwp), intent(in) :: iasta, iaend
  ! local variables
  integer(kind=iwp) :: ia, ja, ka, t, u, v, x, tu, vx, sgn_tu, sgn_vx
  real(kind=wp), allocatable :: f(:)

  call mma_allocate(f,ndeta,label='f')

  do ia=iasta,iaend
    ! f array construction
    f = Zero
    do tu=1,max_ex1a
      t = ex1_a(tu,ia)%p
      u = ex1_a(tu,ia)%q
      sgn_tu = ex1_a(tu,ia)%sgn
      ka = ex1_a(tu,ia)%rank
      f(ka) = f(ka)+sgn_tu*k(t,u)
      do vx=1,max_ex1a
        v = ex1_a(vx,ka)%p
        x = ex1_a(vx,ka)%q
        sgn_vx = ex1_a(vx,ka)%sgn
        ja = ex1_a(vx,ka)%rank
        f(ja) = f(ja)+Half*sgn_tu*sgn_vx*g(v,x,t,u)
      end do
    end do
    ! sigma addition
    ka = 0
    do ja=1,ndeta
      if (f(ja) /= Zero) then
        ka = ka+1
#       ifdef _PROF_
        nflop = nflop+2*ndeta
#       endif
        call daxpy_(ndetb,f(ja),psi(:,ja),1,sgm(:,ia),1)
      end if
    end do
    if (ka > max_ex2a) stop 'exceeded max double excitations'
  end do

  call mma_deallocate(f)

end subroutine sigma2

subroutine sigma3(g,sgm,psi,ibsta,ibend)
! sigma3(ia,ib) = sum_ja,jb sum_tu,vx <jb|E_tu|ib> <ja|E_vx|ia> g_tuvx C(ja,jb)

  ! integrals
  real(kind=wp), intent(in) :: g(my_norb,my_norb,my_norb,my_norb)
  ! wavefunctions
  real(kind=wp), intent(inout) :: sgm(:,:)
  real(kind=wp), intent(in) :: psi(:,:)
  ! determinant indices
  integer(kind=iwp), intent(in) :: ibsta, ibend
  integer(kind=iwp) :: i, n_couples, ib, jb, kb,  &
                       t, u, v, x, & !orbital indices
                       tu, sgn_tu
  integer(kind=iwp), allocatable :: ia(:), ja(:), sgn_vx(:)
  real(kind=wp), allocatable :: f(:), Ctmp(:,:), Vtmp(:)

  call mma_allocate(ja,max_LRs,label='ja')
  call mma_allocate(ia,max_LRs,label='ia')
  call mma_allocate(sgn_vx,max_LRs,label='sgn_vx')
  call mma_allocate(Ctmp,max_LRs,ndetb,label='Ctmp')
  call mma_allocate(Vtmp,max_LRs,label='Vtmp')
  call mma_allocate(f,ndetb,label='f')

  do v=1,my_norb
    do x=1,my_norb
      ! set up L(I), R(I), sgn(I) defined by L(ia) = E_tu R(ia)
      call LRs_init(v,x,nela,my_norb,ja,ia,sgn_vx,n_couples)
      do i=1,n_couples
        Ctmp(i,:) = psi(:,ja(i))*sgn_vx(i)
      end do
      do ib=ibsta,ibend
        f = Zero
        do tu=1,max_ex1b
          t = ex1_b(tu,ib)%p
          u = ex1_b(tu,ib)%q
          sgn_tu = ex1_b(tu,ib)%sgn
          jb = ex1_b(tu,ib)%rank
          ! for singlets, only tu >= vx are used
          if (mult == 1) then
            if (t < v) cycle
            if ((t == v) .and. (u < x)) cycle
            if ((t == v) .and. (u == x)) then
              f(jb) = f(jb)+sgn_tu*Half*g(t,u,v,x)
              cycle
            end if
          end if
          f(jb) = f(jb)+sgn_tu*g(t,u,v,x)
        end do
        ! V(ia) = sum_jb f(jb) C'(ia,jb) for all ia
        Vtmp = Zero
        kb = 0
        ! loop over non-identical excitations
        do tu=1,max_ex1b
          jb = ex1_b(tu,ib)%rank
          if ((jb /= ib) .and. (f(jb) /= Zero)) then
            kb = kb+1
#           ifdef _PROF_
            nflop = nflop+2*n_couples
#           endif
            call daxpy_(n_couples,f(jb),Ctmp(1,jb),1,Vtmp,1)
          end if
        end do
        ! contribution from the identical excitations
        if (f(ib) /= Zero) then
#         ifdef _PROF_
          nflop = nflop+2*n_couples
#         endif
          call daxpy_(n_couples,f(ib),Ctmp(1,ib),1,Vtmp,1)
        end if
        if (kb > max_ex1b) stop 'exceeded max single excitations'
        ! s3(R_ia,ib) = s3(R_ia,ib) + V(ia)
        do i=1,n_couples
          sgm(ia(i),ib) = sgm(ia(i),ib)+Vtmp(i)
        end do
      end do
    end do
  end do

  call mma_deallocate(ja)
  call mma_deallocate(ia)
  call mma_deallocate(sgn_vx)
  call mma_deallocate(Ctmp)
  call mma_deallocate(Vtmp)
  call mma_deallocate(f)

end subroutine sigma3

subroutine ex1_init(k,n,ex1_table)

  use second_quantization, only: binom_coef, ex1, fase, lex_init, lex_next, lexrank

  integer(kind=iwp), intent(in) :: k, n
  type(ex1_struct), intent(out) :: ex1_table(:,:)
  integer(kind=iwp) :: my_ndet, idet, det, tmp, p, q, pq, counter

  !write(u6,'(1x,a)') 'excitation table'
  !write(u6,'(1x,a)') 'p   q   I   J'
  my_ndet = binom_coef(k,n)
  det = lex_init(k,n)
  counter = 0
  do idet=1,my_ndet
    pq = 0
    do p=1,my_norb
      do q=1,my_norb
        tmp = ex1(p,q,det)
        if (tmp /= -1) then
          pq = pq+1
          ex1_table(pq,idet)%p = p
          ex1_table(pq,idet)%q = q
          ex1_table(pq,idet)%sgn = fase(tmp)
          ex1_table(pq,idet)%rank = lexrank(tmp)
          counter = counter+1
          !write(u6,'(1x,4i4)') p,q,idet,fase(tmp)*lexrank(tmp)
        end if
      end do
    end do
    det = lex_next(det)
  end do

end subroutine ex1_init

subroutine LRs_init(p,q,my_nel,my_norb,L,R,sgn,counter)
! for a pair of orbitals p and q, and determinants
! generated by my_nel electrons in my_norb spin orbitals,
! re-enumerate all non-vanishing couples connected
! through: jdet = E_pq idet. The number of couples
! is given by n_det, L(i) = jdet and R(i) = idet,
! where i is a counter from 1 to n_det. The sign of
! jdet is stored in the sgn array.

  use second_quantization, only: binom_coef, ex1, fase, lex_init, lex_next, lexrank

  integer(kind=iwp), intent(in) :: p, q, my_nel, my_norb
  integer(kind=iwp), intent(out) :: L(:), R(:), sgn(:), counter
  integer(kind=iwp) :: my_ndet, idet, det, tmp

  my_ndet = binom_coef(my_nel,my_norb)
  det = lex_init(my_nel,my_norb)
  counter = 0
  do idet=1,my_ndet
    tmp = ex1(p,q,det)
    if (tmp /= -1) then
      counter = counter+1
      R(counter) = idet
      sgn(counter) = fase(tmp)
      L(counter) = lexrank(tmp)
    end if
    det = lex_next(det)
  end do

end subroutine LRs_init

! Extensions to mma_interfaces, using preprocessor templates
! (see src/mma_util/stdalloc.f)

! Define ex1_cptr2loff, ex1_mma_allo_2D, ex1_mma_allo_2D_lim, ex1_mma_free_2D
#define _TYPE_ type(ex1_struct)
#  define _FUNC_NAME_ ex1_cptr2loff
#  include "cptr2loff_template.fh"
#  undef _FUNC_NAME_
#  define _SUBR_NAME_ ex1_mma
#  define _DIMENSIONS_ 2
#  define _DEF_LABEL_ 'ex1_mma'
#  include "mma_allo_template.fh"
#  undef _SUBR_NAME_
#  undef _DIMENSIONS_
#  undef _DEF_LABEL_
#undef _TYPE_

end module faroald
