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
use Constants, only: Zero, One, Half
use Definitions, only: wp, iwp, u6
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

! integral storage in Faroald format
real(kind=wp), allocatable:: gtuvx(:,:,:,:), htu(:,:)
#ifdef _PROF_
integer(kind=int64) :: nflop
#endif

public :: ex1_a, ex1_b, ex1_init, max_LRs, max_ex1a, max_ex1b, max_ex2a, max_ex2b, mult, my_ndet, my_nel, my_norb, ndeta, ndetb, &
          nela, nelb, nhoa, nhob, sigma_update, gtuvx, htu, one_pdm, transition_one_pdm, two_pdm, transition_two_pdm, fold_two_pdm

! Extensions to mma interfaces

interface mma_allocate
  module procedure :: ex1_mma_allo_2D, ex1_mma_allo_2D_lim
end interface
interface mma_deallocate
  module procedure :: ex1_mma_free_2D
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
  call gadgop(sgm,ndeta*ndetb,'+')

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
    write(u6,'(1x,a,2(f10.3,a))') 'sigma update: ',walltime,' s, ',flops*1.0e-9_wp,' Gflops.'
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
        sgm(1:ndeta,ib) = sgm(1:ndeta,ib)+f(jb)*psi(1:ndeta,jb)
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
        sgm(1:ndetb,ia) = sgm(1:ndetb,ia)+f(ja)*psi(:,ja)
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
  integer(kind=iwp) :: i, n_couples, ib, jb, kb, &
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
            Vtmp(1:n_couples) = Vtmp(1:n_couples)+f(jb)*Ctmp(1:n_couples,jb)
          end if
        end do
        ! contribution from the identical excitations
        if (f(ib) /= Zero) then
#         ifdef _PROF_
          nflop = nflop+2*n_couples
#         endif
          Vtmp(1:n_couples) = Vtmp(1:n_couples)+f(ib)*Ctmp(1:n_couples,ib)
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

subroutine one_pdm(psi, d, sd, da, db)
! Compute the spin-summed one-particle density matrix and spin-density
! matrix for a CI vector in determinant product representation.
!
!   da(p,q) = <Psi| a^+_{p,alpha} a_{q,alpha} |Psi>
!   db(p,q) = <Psi| a^+_{p,beta } a_{q,beta } |Psi>
!
!   d (p,q) = da(p,q) + db(p,q)
!   sd(p,q) = da(p,q) - db(p,q)
!
! If the CI vector is not normalized, the densities are scaled by
! <Psi|Psi>.

  real(kind=wp), intent(in)  :: psi(:,:)
  real(kind=wp), intent(out) :: d(my_norb,my_norb)
  real(kind=wp), intent(out) :: sd(my_norb,my_norb)
  real(kind=wp), intent(out), optional :: da(my_norb,my_norb)
  real(kind=wp), intent(out), optional :: db(my_norb,my_norb)

  call transition_one_pdm(psi,psi,d,sd,da,db)

end subroutine one_pdm


subroutine transition_one_pdm(bra, ket, d, sd, da, db)
! Compute the spin-resolved transition one-particle density matrices
! between two CI vectors:
!
!   da(p,q) = <bra| a^+_{p,alpha} a_{q,alpha} |ket>
!   db(p,q) = <bra| a^+_{p,beta } a_{q,beta } |ket>
!
! and return
!
!   d (p,q) = da(p,q) + db(p,q)
!   sd(p,q) = da(p,q) - db(p,q)
!
! For an ordinary density matrix, call this routine with bra == ket.
!
! The determinant excitation information is taken from ex1_a and ex1_b.
! Each table entry represents
!
!   |I> = E_pq |J>
!
! with determinant rank I and fermionic sign sgn.

  real(kind=wp), intent(in)  :: bra(:,:)
  real(kind=wp), intent(in)  :: ket(:,:)
  real(kind=wp), intent(out) :: d(my_norb,my_norb)
  real(kind=wp), intent(out) :: sd(my_norb,my_norb)
  real(kind=wp), intent(out), optional :: da(my_norb,my_norb)
  real(kind=wp), intent(out), optional :: db(my_norb,my_norb)

  real(kind=wp), allocatable :: da_loc(:,:), db_loc(:,:)

  call mma_allocate(da_loc,my_norb,my_norb,label='da_loc')
  call mma_allocate(db_loc,my_norb,my_norb,label='db_loc')

  call alpha_transition_one_pdm(bra,ket,da_loc)
  call beta_transition_one_pdm (bra,ket,db_loc)

  d  = da_loc + db_loc
  sd = da_loc - db_loc

  if (present(da)) da = da_loc
  if (present(db)) db = db_loc

  call mma_deallocate(da_loc)
  call mma_deallocate(db_loc)

end subroutine transition_one_pdm

subroutine alpha_transition_one_pdm(bra, ket, da)
! Alpha-spin contribution to the transition one-particle density matrix:
!
!   da(p,q) = sum_ia,ja,ib bra(ia,ib)
!             <ia|E_pq|ja> ket(ja,ib)

  real(kind=wp), intent(in)  :: bra(:,:)
  real(kind=wp), intent(in)  :: ket(:,:)
  real(kind=wp), intent(out) :: da(my_norb,my_norb)

  integer(kind=iwp) :: ja, ia, jasta, jaend
  integer(kind=iwp) :: pq, p, q, sgn
  real(kind=wp) :: contribution

  da = Zero

  call par_range(ndeta,jasta,jaend)

  do ja = jasta, jaend
    do pq = 1, max_ex1a

      p   = ex1_a(pq,ja)%p
      q   = ex1_a(pq,ja)%q
      sgn = ex1_a(pq,ja)%sgn
      ia  = ex1_a(pq,ja)%rank

      contribution = dot_product(bra(ia,1:ndetb),ket(ja,1:ndetb))

      da(p,q) = da(p,q) + real(sgn,kind=wp)*contribution

    end do
  end do

  call gadgop(da,my_norb*my_norb,'+')

end subroutine alpha_transition_one_pdm

subroutine beta_transition_one_pdm(bra, ket, db)
! Beta-spin contribution to the transition one-particle density matrix:
!
!   db(p,q) = sum_ib,jb,ia bra(ia,ib)
!             <ib|E_pq|jb> ket(ia,jb)
!
! The excitation table ex1_b stores, for each source determinant jb,
! the destination determinant ib, orbital pair p,q, and phase.

  real(kind=wp), intent(in)  :: bra(:,:)
  real(kind=wp), intent(in)  :: ket(:,:)
  real(kind=wp), intent(out) :: db(my_norb,my_norb)

  integer(kind=iwp) :: jb, ib, jbsta, jbend
  integer(kind=iwp) :: pq, p, q, sgn
  real(kind=wp) :: contribution

  db = Zero

  ! In analogy with sigma_update, distribute the determinant work if
  ! running under the existing parallel environment.
  call par_range(ndetb,jbsta,jbend)

  do jb = jbsta, jbend
    do pq = 1, max_ex1b

      p   = ex1_b(pq,jb)%p
      q   = ex1_b(pq,jb)%q
      sgn = ex1_b(pq,jb)%sgn
      ib  = ex1_b(pq,jb)%rank

      contribution = dot_product(bra(1:ndeta,ib),ket(1:ndeta,jb))

      db(p,q) = db(p,q) + real(sgn,kind=wp)*contribution

    end do
  end do

  call gadgop(db,my_norb*my_norb,'+')

end subroutine beta_transition_one_pdm

subroutine two_pdm(psi, p2, p2prod)
! Compute the spin-free two-particle density matrix for one CI vector.
!
! The returned p2 is the normal-ordered spin-free two-particle density
!
!   p2(t,u,v,x) = <Psi| E_tu E_vx - delta(u,v) E_tx |Psi>
!
! where
!
!   E_tu = sum_sigma a^+_{t sigma} a_{u sigma}.
!
! If p2prod is present, it receives the unnormal-ordered product density
!
!   p2prod(t,u,v,x) = <Psi| E_tu E_vx |Psi>.
!
! If psi is not normalized, the density matrices are scaled by <Psi|Psi>.

  real(kind=wp), intent(in)  :: psi(:,:)
  real(kind=wp), intent(out) :: p2(my_norb,my_norb,my_norb,my_norb)
  real(kind=wp), intent(out), optional :: p2prod(my_norb,my_norb,my_norb,my_norb)

  call transition_two_pdm(psi,psi,p2,p2prod)

end subroutine two_pdm


subroutine transition_two_pdm(bra, ket, p2, p2prod)
! Compute the spin-free transition two-particle density matrix.
!
! The returned p2 is
!
!   p2(t,u,v,x) = <bra| E_tu E_vx - delta(u,v) E_tx |ket>.
!
! If p2prod is present, it receives
!
!   p2prod(t,u,v,x) = <bra| E_tu E_vx |ket>.
!
! The construction first forms the product density <E_tu E_vx> by
! applying two one-particle excitation operators through the existing
! alpha and beta excitation tables. Then the contraction
!
!   delta(u,v) <E_tx>
!
! is subtracted using transition_one_pdm.

  real(kind=wp), intent(in)  :: bra(:,:)
  real(kind=wp), intent(in)  :: ket(:,:)
  real(kind=wp), intent(out) :: p2(my_norb,my_norb,my_norb,my_norb)
  real(kind=wp), intent(out), optional :: p2prod(my_norb,my_norb,my_norb,my_norb)

  real(kind=wp), allocatable :: pprod(:,:,:,:)
  real(kind=wp), allocatable :: d1(:,:), sd1(:,:)
  integer(kind=iwp) :: t, u, v, x

  call mma_allocate(pprod,my_norb,my_norb,my_norb,my_norb,label='pprod')
  call mma_allocate(d1,my_norb,my_norb,label='d1')
  call mma_allocate(sd1,my_norb,my_norb,label='sd1')

  call transition_two_pdm_product(bra,ket,pprod)

  ! One-particle transition density needed for the contraction term.
  call transition_one_pdm(bra,ket,d1,sd1)

  p2 = pprod

  do x = 1, my_norb
    do v = 1, my_norb
      do u = 1, my_norb
        do t = 1, my_norb
          if (u == v) p2(t,u,v,x) = p2(t,u,v,x) - d1(t,x)
        end do
      end do
    end do
  end do

  if (present(p2prod)) p2prod = pprod

  call mma_deallocate(sd1)
  call mma_deallocate(d1)
  call mma_deallocate(pprod)

end subroutine transition_two_pdm


subroutine transition_two_pdm_product(bra, ket, pprod)
! Compute the spin-free product density
!
!   pprod(t,u,v,x) = <bra| E_tu E_vx |ket>
!
! without subtracting the contraction delta(u,v) <E_tx>.
!
! The spin-free product is built from four spin blocks:
!
!   alpha-alpha : E_tu(alpha) E_vx(alpha)
!   beta-beta   : E_tu(beta ) E_vx(beta )
!   alpha-beta  : E_tu(alpha) E_vx(beta )
!   beta-alpha  : E_tu(beta ) E_vx(alpha)
!
! The rightmost operator E_vx acts first.

  real(kind=wp), intent(in)  :: bra(:,:)
  real(kind=wp), intent(in)  :: ket(:,:)
  real(kind=wp), intent(out) :: pprod(my_norb,my_norb,my_norb,my_norb)

  pprod = Zero

  call alpha_alpha_two_pdm_product(bra,ket,pprod)
  call beta_beta_two_pdm_product  (bra,ket,pprod)
  call alpha_beta_two_pdm_product (bra,ket,pprod)
  call beta_alpha_two_pdm_product (bra,ket,pprod)

end subroutine transition_two_pdm_product


subroutine alpha_alpha_two_pdm_product(bra, ket, pprod)
! Add the alpha-alpha contribution
!
!   <bra| E_tu(alpha) E_vx(alpha) |ket>
!
! to pprod(t,u,v,x). The rightmost operator E_vx(alpha) acts first:
!
!   ja --E_vx--> ka --E_tu--> ia

  real(kind=wp), intent(in)    :: bra(:,:)
  real(kind=wp), intent(in)    :: ket(:,:)
  real(kind=wp), intent(inout) :: pprod(my_norb,my_norb,my_norb,my_norb)

  integer(kind=iwp) :: ja, ka, ia
  integer(kind=iwp) :: tu, vx
  integer(kind=iwp) :: t, u, v, x
  integer(kind=iwp) :: sgn_tu, sgn_vx
  real(kind=wp) :: contribution

  do ja = 1, ndeta
    do vx = 1, max_ex1a

      v      = ex1_a(vx,ja)%p
      x      = ex1_a(vx,ja)%q
      sgn_vx = ex1_a(vx,ja)%sgn
      ka     = ex1_a(vx,ja)%rank

      do tu = 1, max_ex1a

        t      = ex1_a(tu,ka)%p
        u      = ex1_a(tu,ka)%q
        sgn_tu = ex1_a(tu,ka)%sgn
        ia     = ex1_a(tu,ka)%rank

        contribution = dot_product(bra(ia,1:ndetb),ket(ja,1:ndetb))

        pprod(t,u,v,x) = pprod(t,u,v,x)                                 &
             + real(sgn_tu*sgn_vx,kind=wp)*contribution

      end do
    end do
  end do

end subroutine alpha_alpha_two_pdm_product


subroutine beta_beta_two_pdm_product(bra, ket, pprod)
! Add the beta-beta contribution
!
!   <bra| E_tu(beta) E_vx(beta) |ket>
!
! to pprod(t,u,v,x). The rightmost operator E_vx(beta) acts first:
!
!   jb --E_vx--> kb --E_tu--> ib

  real(kind=wp), intent(in)    :: bra(:,:)
  real(kind=wp), intent(in)    :: ket(:,:)
  real(kind=wp), intent(inout) :: pprod(my_norb,my_norb,my_norb,my_norb)

  integer(kind=iwp) :: jb, kb, ib
  integer(kind=iwp) :: tu, vx
  integer(kind=iwp) :: t, u, v, x
  integer(kind=iwp) :: sgn_tu, sgn_vx
  real(kind=wp) :: contribution

  do jb = 1, ndetb
    do vx = 1, max_ex1b

      v      = ex1_b(vx,jb)%p
      x      = ex1_b(vx,jb)%q
      sgn_vx = ex1_b(vx,jb)%sgn
      kb     = ex1_b(vx,jb)%rank

      do tu = 1, max_ex1b

        t      = ex1_b(tu,kb)%p
        u      = ex1_b(tu,kb)%q
        sgn_tu = ex1_b(tu,kb)%sgn
        ib     = ex1_b(tu,kb)%rank

        contribution = dot_product(bra(1:ndeta,ib),ket(1:ndeta,jb))

        pprod(t,u,v,x) = pprod(t,u,v,x)                                 &
             + real(sgn_tu*sgn_vx,kind=wp)*contribution

      end do
    end do
  end do

end subroutine beta_beta_two_pdm_product


subroutine alpha_beta_two_pdm_product(bra, ket, pprod)
! Add the alpha-beta contribution
!
!   <bra| E_tu(alpha) E_vx(beta) |ket>
!
! to pprod(t,u,v,x). Since alpha and beta strings are stored separately,
! the excitation signs factor into the alpha and beta signs:
!
!   ja --E_tu(alpha)--> ia
!   jb --E_vx(beta )--> ib

  real(kind=wp), intent(in)    :: bra(:,:)
  real(kind=wp), intent(in)    :: ket(:,:)
  real(kind=wp), intent(inout) :: pprod(my_norb,my_norb,my_norb,my_norb)

  integer(kind=iwp) :: ja, ia, jb, ib
  integer(kind=iwp) :: tu, vx
  integer(kind=iwp) :: t, u, v, x
  integer(kind=iwp) :: sgn_tu, sgn_vx

  do ja = 1, ndeta
    do tu = 1, max_ex1a

      t      = ex1_a(tu,ja)%p
      u      = ex1_a(tu,ja)%q
      sgn_tu = ex1_a(tu,ja)%sgn
      ia     = ex1_a(tu,ja)%rank

      do jb = 1, ndetb
        do vx = 1, max_ex1b

          v      = ex1_b(vx,jb)%p
          x      = ex1_b(vx,jb)%q
          sgn_vx = ex1_b(vx,jb)%sgn
          ib     = ex1_b(vx,jb)%rank

          pprod(t,u,v,x) = pprod(t,u,v,x)                               &
               + real(sgn_tu*sgn_vx,kind=wp)*bra(ia,ib)*ket(ja,jb)

        end do
      end do

    end do
  end do

end subroutine alpha_beta_two_pdm_product


subroutine beta_alpha_two_pdm_product(bra, ket, pprod)
! Add the beta-alpha contribution
!
!   <bra| E_tu(beta) E_vx(alpha) |ket>
!
! to pprod(t,u,v,x):
!
!   jb --E_tu(beta )--> ib
!   ja --E_vx(alpha)--> ia

  real(kind=wp), intent(in)    :: bra(:,:)
  real(kind=wp), intent(in)    :: ket(:,:)
  real(kind=wp), intent(inout) :: pprod(my_norb,my_norb,my_norb,my_norb)

  integer(kind=iwp) :: ja, ia, jb, ib
  integer(kind=iwp) :: tu, vx
  integer(kind=iwp) :: t, u, v, x
  integer(kind=iwp) :: sgn_tu, sgn_vx

  do jb = 1, ndetb
    do tu = 1, max_ex1b

      t      = ex1_b(tu,jb)%p
      u      = ex1_b(tu,jb)%q
      sgn_tu = ex1_b(tu,jb)%sgn
      ib     = ex1_b(tu,jb)%rank

      do ja = 1, ndeta
        do vx = 1, max_ex1a

          v      = ex1_a(vx,ja)%p
          x      = ex1_a(vx,ja)%q
          sgn_vx = ex1_a(vx,ja)%sgn
          ia     = ex1_a(vx,ja)%rank

          pprod(t,u,v,x) = pprod(t,u,v,x)                               &
               + real(sgn_tu*sgn_vx,kind=wp)*bra(ia,ib)*ket(ja,jb)

        end do
      end do

    end do
  end do

end subroutine beta_alpha_two_pdm_product

subroutine fold_two_pdm(p2, p2_fold, average)
! Fold a full four-index spin-free two-particle density matrix
!
!   p2(t,u,v,x)
!
! into packed pair-pair triangular storage:
!
!   p2_fold(ij)
!
! where
!
!   tu = pair_index(t,u),  t >= u
!   vx = pair_index(v,x),  v >= x
!   ij = pair_index(tu,vx), tu >= vx
!
! Thus p2_fold has length
!
!   npair2 = npair*(npair+1)/2
!   npair  = my_norb*(my_norb+1)/2
!
! By default, this routine performs a direct canonical fold:
!
!   p2_fold((tu,vx)) = p2(t,u,v,x)
!
! using only canonical representatives t>=u, v>=x, tu>=vx.
!
! If average=.true., all tensor elements mapping to the same packed
! position are averaged. This is useful as a diagnostic if the full
! tensor has small numerical deviations from the expected packed
! symmetries.

  real(kind=wp), intent(in)  :: p2(my_norb,my_norb,my_norb,my_norb)
  real(kind=wp), intent(out) :: p2_fold(:)
  logical, intent(in), optional :: average

  integer(kind=iwp) :: t, u, v, x
  integer(kind=iwp) :: tu, vx, tuvx
  integer(kind=iwp) :: npair, npair2
  logical :: do_average
  real(kind=wp), allocatable :: weight(:)

  npair  = my_norb*(my_norb+1)/2
  npair2 = npair*(npair+1)/2

  if (size(p2_fold) < npair2) then
    Write (u6,*) 'fold_two_pdm: p2_fold too small'
    Call Abend()
  end if

  do_average = .false.
  if (present(average)) do_average = average

  p2_fold(:) = Zero

  if (.not. do_average) then

    ! Canonical Fold2-like packing.
    !
    ! Only canonical orbital-pair representatives are used:
    !
    !   t >= u
    !   v >= x
    !   pair(t,u) >= pair(v,x)

    do t = 1, my_norb
      do u = 1, t

        tu = faroald_pair_index(t,u)

        do v = 1, my_norb
          do x = 1, v

            vx = faroald_pair_index(v,x)

            if (tu < vx) cycle

            tuvx = faroald_pair_index(tu,vx)

            p2_fold(tuvx) = p2(t,u,v,x)

          end do
        end do

      end do
    end do

  else

    ! Symmetry-averaged packing.
    !
    ! Every full tensor element is mapped to the same packed pair-pair
    ! address as its symmetry-equivalent partners. The stored value is
    ! the arithmetic average over all entries mapped to that address.

    call mma_allocate(weight,npair2,label='fold_two_pdm_weight')
    weight(:) = Zero

    do t = 1, my_norb
      do u = 1, my_norb

        tu = faroald_pair_index(t,u)

        do v = 1, my_norb
          do x = 1, my_norb

            vx = faroald_pair_index(v,x)
            tuvx = faroald_pair_index(tu,vx)

            p2_fold(tuvx) = p2_fold(tuvx) + p2(t,u,v,x)
            weight(tuvx)  = weight(tuvx)  + One

          end do
        end do

      end do
    end do

    do tuvx = 1, npair2
      if (weight(tuvx) /= Zero) p2_fold(tuvx) = p2_fold(tuvx)/weight(tuvx)
    end do

    call mma_deallocate(weight)

  end if

end subroutine fold_two_pdm


integer(kind=iwp) function faroald_pair_index(i, j) result(ij)
! Return the packed lower-triangular index for an unordered pair (i,j).
!
! If i >= j:
!
!   ij = i*(i-1)/2 + j
!
! Otherwise the pair is swapped.

  integer(kind=iwp), intent(in) :: i, j

  if (i >= j) then
    ij = i*(i-1)/2 + j
  else
    ij = j*(j-1)/2 + i
  end if

end function faroald_pair_index

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
! (see mma_util/stdalloc.F90)

! Define ex1_mma_allo_2D, ex1_mma_allo_2D_lim, ex1_mma_free_2D
#define _TYPE_ type(ex1_struct)
#  define _SUBR_NAME_ ex1_mma
#  define _DIMENSIONS_ 2
#  define _DEF_LABEL_ 'ex1_mma'
#  include "mma_allo_template.fh"
#  undef _SUBR_NAME_
#  undef _DIMENSIONS_
#  undef _DEF_LABEL_
#undef _TYPE_

end module faroald
