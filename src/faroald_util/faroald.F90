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

use second_quantization, only: binom_coef
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, Half
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
real(kind=wp), allocatable :: gtuvx(:,:,:,:), htu(:,:)

#ifdef _PROF_
integer(kind=int64) :: nflop
#endif

integer(kind=iwp) :: npat

integer(kind=iwp), allocatable :: occpat(:,:),  &    ! (my_norb,npat)
                                  ipat_of_det(:)
integer(kind=iwp), allocatable :: ndoub(:), nsing(:)
!
! nspin_comb : raw alpha/beta assignments among open shells
! ncomb      : Lucia-compatible counting of diagonal states
!
integer(kind=iwp), allocatable :: nspin_comb(:), nComb(:)
integer(kind=iwp), allocatable :: ibcomb(:)
integer(kind=iwp) :: ncomb_tot
integer(kind=iwp), allocatable :: ictsdt(:)

public :: ex1_a, ex1_b, ex1_init, fold_two_pdm, gtuvx, htu, max_ex1a, max_ex1b, max_ex2a, max_ex2b, max_LRs, mult, my_ndet, &
          my_nel, my_norb, ndeta, ndetb, nela, nelb, nhoa, nhob, one_pdm, sigma_update, transition_one_pdm, transition_two_pdm, &
          two_pdm, hDiag, verify_occ_patterns, build_patterns, analyse_patterns, ndoub, nsing, build_pattern_diagonal, &
          build_pattern_diagonal_approx, nSpin_Comb, nComb, ibComb, nComb_tot, pattern_combinations
public :: npat, occpat, ipat_of_det
public :: combination_occupations
public :: combination_diagonal
public :: build_ictsdt, ictsdt
public :: faroald_epq_psi, traci_faroald

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

subroutine one_pdm(psi,d,sd,da,db)
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

  real(kind=wp), intent(in) :: psi(:,:)
  real(kind=wp), intent(out) :: d(my_norb,my_norb), sd(my_norb,my_norb)
  real(kind=wp), intent(out), optional :: da(my_norb,my_norb), db(my_norb,my_norb)

  call transition_one_pdm(psi,psi,d,sd,da,db)

end subroutine one_pdm

subroutine transition_one_pdm(bra,ket,d,sd,da,db)
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

  real(kind=wp), intent(in) :: bra(:,:), ket(:,:)
  real(kind=wp), intent(out) :: d(my_norb,my_norb), sd(my_norb,my_norb)
  real(kind=wp), intent(out), optional :: da(my_norb,my_norb), db(my_norb,my_norb)
  real(kind=wp), allocatable :: da_loc(:,:), db_loc(:,:)

  call mma_allocate(da_loc,my_norb,my_norb,label='da_loc')
  call mma_allocate(db_loc,my_norb,my_norb,label='db_loc')

  call alpha_transition_one_pdm(bra,ket,da_loc)
  call beta_transition_one_pdm(bra,ket,db_loc)

  d(:,:) = da_loc(:,:)+db_loc(:,:)
  sd(:,:) = da_loc(:,:)-db_loc(:,:)

  if (present(da)) da(:,:) = da_loc(:,:)
  if (present(db)) db(:,:) = db_loc(:,:)

  call mma_deallocate(da_loc)
  call mma_deallocate(db_loc)

end subroutine transition_one_pdm

subroutine alpha_transition_one_pdm(bra,ket,da)
! Alpha-spin contribution to the transition one-particle density matrix:
!
!   da(p,q) = sum_ia,ja,ib bra(ia,ib)
!             <ia|E_pq|ja> ket(ja,ib)

  real(kind=wp), intent(in) :: bra(:,:), ket(:,:)
  real(kind=wp), intent(out) :: da(my_norb,my_norb)
  integer(kind=iwp) :: ia, ja, jaend, jasta, p, pq, q, sgn
  real(kind=wp) :: contribution

  da = Zero

  call par_range(ndeta,jasta,jaend)

  do ja=jasta,jaend
    do pq=1,max_ex1a

      p = ex1_a(pq,ja)%p
      q = ex1_a(pq,ja)%q
      sgn = ex1_a(pq,ja)%sgn
      ia = ex1_a(pq,ja)%rank

      contribution = dot_product(bra(ia,1:ndetb),ket(ja,1:ndetb))

      da(p,q) = da(p,q)+real(sgn,kind=wp)*contribution

    end do
  end do

  call gadgop(da,my_norb*my_norb,'+')

end subroutine alpha_transition_one_pdm

subroutine beta_transition_one_pdm(bra,ket,db)
! Beta-spin contribution to the transition one-particle density matrix:
!
!   db(p,q) = sum_ib,jb,ia bra(ia,ib)
!             <ib|E_pq|jb> ket(ia,jb)
!
! The excitation table ex1_b stores, for each source determinant jb,
! the destination determinant ib, orbital pair p,q, and phase.

  real(kind=wp), intent(in) :: bra(:,:), ket(:,:)
  real(kind=wp), intent(out) :: db(my_norb,my_norb)
  integer(kind=iwp) :: ib, jb, jbend, jbsta, p, pq, q, sgn
  real(kind=wp) :: contribution

  db = Zero

  ! In analogy with sigma_update, distribute the determinant work if
  ! running under the existing parallel environment.
  call par_range(ndetb,jbsta,jbend)

  do jb=jbsta,jbend
    do pq=1,max_ex1b

      p = ex1_b(pq,jb)%p
      q = ex1_b(pq,jb)%q
      sgn = ex1_b(pq,jb)%sgn
      ib = ex1_b(pq,jb)%rank

      contribution = dot_product(bra(1:ndeta,ib),ket(1:ndeta,jb))

      db(p,q) = db(p,q)+real(sgn,kind=wp)*contribution

    end do
  end do

  call gadgop(db,my_norb*my_norb,'+')

end subroutine beta_transition_one_pdm

subroutine two_pdm(psi,p2,p2prod)
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

  real(kind=wp), intent(in) :: psi(:,:)
  real(kind=wp), intent(out) :: p2(my_norb,my_norb,my_norb,my_norb)
  real(kind=wp), intent(out), optional :: p2prod(my_norb,my_norb,my_norb,my_norb)

  call transition_two_pdm(psi,psi,p2,p2prod)

end subroutine two_pdm

subroutine transition_two_pdm(bra,ket,p2,p2prod)
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

  real(kind=wp), intent(in) :: bra(:,:), ket(:,:)
  real(kind=wp), intent(out) :: p2(my_norb,my_norb,my_norb,my_norb)
  real(kind=wp), intent(out), optional :: p2prod(my_norb,my_norb,my_norb,my_norb)
  real(kind=wp), allocatable :: d1(:,:), pprod(:,:,:,:), sd1(:,:)
  integer(kind=iwp) :: u

  call mma_allocate(pprod,my_norb,my_norb,my_norb,my_norb,label='pprod')
  call mma_allocate(d1,my_norb,my_norb,label='d1')
  call mma_allocate(sd1,my_norb,my_norb,label='sd1')

  call transition_two_pdm_product(bra,ket,pprod)

  ! One-particle transition density needed for the contraction term.
  call transition_one_pdm(bra,ket,d1,sd1)

  p2(:,:,:,:) = pprod(:,:,:,:)

  do u=1,my_norb
    p2(:,u,u,:) = p2(:,u,u,:)-d1(:,:)
  end do

  if (present(p2prod)) p2prod(:,:,:,:) = pprod(:,:,:,:)

  call mma_deallocate(sd1)
  call mma_deallocate(d1)
  call mma_deallocate(pprod)

end subroutine transition_two_pdm

subroutine transition_two_pdm_product(bra,ket,pprod)
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

  real(kind=wp), intent(in) :: bra(:,:), ket(:,:)
  real(kind=wp), intent(out) :: pprod(my_norb,my_norb,my_norb,my_norb)

  pprod(:,:,:,:) = Zero

  call alpha_alpha_two_pdm_product(bra,ket,pprod)
  call beta_beta_two_pdm_product(bra,ket,pprod)
  call alpha_beta_two_pdm_product(bra,ket,pprod)
  call beta_alpha_two_pdm_product(bra,ket,pprod)

end subroutine transition_two_pdm_product

subroutine alpha_alpha_two_pdm_product(bra,ket,pprod)
! Add the alpha-alpha contribution
!
!   <bra| E_tu(alpha) E_vx(alpha) |ket>
!
! to pprod(t,u,v,x). The rightmost operator E_vx(alpha) acts first:
!
!   ja --E_vx--> ka --E_tu--> ia

  real(kind=wp), intent(in) :: bra(:,:), ket(:,:)
  real(kind=wp), intent(inout) :: pprod(my_norb,my_norb,my_norb,my_norb)
  integer(kind=iwp) :: ia, ja, ka, sgn_tu, sgn_vx, t, tu, u, v, vx, x
  real(kind=wp) :: contribution

  do ja=1,ndeta
    do vx=1,max_ex1a

      v = ex1_a(vx,ja)%p
      x = ex1_a(vx,ja)%q
      sgn_vx = ex1_a(vx,ja)%sgn
      ka = ex1_a(vx,ja)%rank

      do tu=1,max_ex1a

        t = ex1_a(tu,ka)%p
        u = ex1_a(tu,ka)%q
        sgn_tu = ex1_a(tu,ka)%sgn
        ia = ex1_a(tu,ka)%rank

        contribution = dot_product(bra(ia,1:ndetb),ket(ja,1:ndetb))

        pprod(t,u,v,x) = pprod(t,u,v,x)+real(sgn_tu*sgn_vx,kind=wp)*contribution

      end do
    end do
  end do

end subroutine alpha_alpha_two_pdm_product

subroutine beta_beta_two_pdm_product(bra,ket,pprod)
! Add the beta-beta contribution
!
!   <bra| E_tu(beta) E_vx(beta) |ket>
!
! to pprod(t,u,v,x). The rightmost operator E_vx(beta) acts first:
!
!   jb --E_vx--> kb --E_tu--> ib

  real(kind=wp), intent(in) :: bra(:,:), ket(:,:)
  real(kind=wp), intent(inout) :: pprod(my_norb,my_norb,my_norb,my_norb)
  integer(kind=iwp) :: ib, jb, kb, sgn_tu, sgn_vx, t, tu, u, v, vx, x
  real(kind=wp) :: contribution

  do jb=1,ndetb
    do vx=1,max_ex1b

      v = ex1_b(vx,jb)%p
      x = ex1_b(vx,jb)%q
      sgn_vx = ex1_b(vx,jb)%sgn
      kb = ex1_b(vx,jb)%rank

      do tu=1,max_ex1b

        t = ex1_b(tu,kb)%p
        u = ex1_b(tu,kb)%q
        sgn_tu = ex1_b(tu,kb)%sgn
        ib = ex1_b(tu,kb)%rank

        contribution = dot_product(bra(1:ndeta,ib),ket(1:ndeta,jb))

        pprod(t,u,v,x) = pprod(t,u,v,x)+real(sgn_tu*sgn_vx,kind=wp)*contribution

      end do
    end do
  end do

end subroutine beta_beta_two_pdm_product

subroutine alpha_beta_two_pdm_product(bra,ket,pprod)
! Add the alpha-beta contribution
!
!   <bra| E_tu(alpha) E_vx(beta) |ket>
!
! to pprod(t,u,v,x). Since alpha and beta strings are stored separately,
! the excitation signs factor into the alpha and beta signs:
!
!   ja --E_tu(alpha)--> ia
!   jb --E_vx(beta )--> ib

  real(kind=wp), intent(in) :: bra(:,:), ket(:,:)
  real(kind=wp), intent(inout) :: pprod(my_norb,my_norb,my_norb,my_norb)
  integer(kind=iwp) :: ia, ib, ja, jb, sgn_tu, sgn_vx, t, tu, u, v, vx, x

  do ja=1,ndeta
    do tu=1,max_ex1a

      t = ex1_a(tu,ja)%p
      u = ex1_a(tu,ja)%q
      sgn_tu = ex1_a(tu,ja)%sgn
      ia = ex1_a(tu,ja)%rank

      do jb=1,ndetb
        do vx=1,max_ex1b

          v = ex1_b(vx,jb)%p
          x = ex1_b(vx,jb)%q
          sgn_vx = ex1_b(vx,jb)%sgn
          ib = ex1_b(vx,jb)%rank

          pprod(t,u,v,x) = pprod(t,u,v,x)+real(sgn_tu*sgn_vx,kind=wp)*bra(ia,ib)*ket(ja,jb)

        end do
      end do

    end do
  end do

end subroutine alpha_beta_two_pdm_product

subroutine beta_alpha_two_pdm_product(bra,ket,pprod)
! Add the beta-alpha contribution
!
!   <bra| E_tu(beta) E_vx(alpha) |ket>
!
! to pprod(t,u,v,x):
!
!   jb --E_tu(beta )--> ib
!   ja --E_vx(alpha)--> ia

  real(kind=wp), intent(in) :: bra(:,:), ket(:,:)
  real(kind=wp), intent(inout) :: pprod(my_norb,my_norb,my_norb,my_norb)
  integer(kind=iwp) :: ia, ib, ja, jb, sgn_tu, sgn_vx, t, tu, u, v, vx, x

  do jb=1,ndetb
    do tu=1,max_ex1b

      t = ex1_b(tu,jb)%p
      u = ex1_b(tu,jb)%q
      sgn_tu = ex1_b(tu,jb)%sgn
      ib = ex1_b(tu,jb)%rank

      do ja=1,ndeta
        do vx=1,max_ex1a

          v = ex1_a(vx,ja)%p
          x = ex1_a(vx,ja)%q
          sgn_vx = ex1_a(vx,ja)%sgn
          ia = ex1_a(vx,ja)%rank

          pprod(t,u,v,x) = pprod(t,u,v,x)+real(sgn_tu*sgn_vx,kind=wp)*bra(ia,ib)*ket(ja,jb)

        end do
      end do

    end do
  end do

end subroutine beta_alpha_two_pdm_product

subroutine fold_two_pdm(p2,p2_fold,p2a_fold)
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

  use Index_Functions, only: iTri, nTri_Elem
  use Constants, only: Half

  real(kind=wp), intent(in) :: p2(my_norb,my_norb,my_norb,my_norb)
  real(kind=wp), intent(out) :: p2_fold(:), p2a_fold(:)
  integer(kind=iwp) :: npair, npair2, t, tu, tuvx, u, v, vx, x, x_max

  npair = nTri_Elem(my_norb)
  npair2 = nTri_Elem(npair)

  if (size(p2_fold) < npair2) then
    write(u6,*) 'fold_two_pdm: p2_fold too small'
    call Abend()
  end if

  p2_fold(:) = Zero
  p2a_fold(:) = Zero

  ! Canonical Fold2-like packing.
  !
  ! Only canonical orbital-pair representatives are used:
  !
  !   t >= u
  !   v >= x
  !   pair(t,u) >= pair(v,x)

  do t=1,my_norb
    do u=1,t

      tu = iTri(t,u)

      do v=1,t
        x_max = v
        if (v == t) x_max = u
        do x=1,x_max

          vx = iTri(v,x)

          tuvx = iTri(tu,vx)

          if (v == x) then
            p2_fold(tuvx) = p2(t,u,v,x)
          else
            p2_fold(tuvx) = p2(t,u,v,x)+p2(t,u,x,v)
          end if

          if ((t /= u) .and. (v /= x)) p2a_fold(tuvx) = p2(t,u,v,x)-p2(t,u,x,v)

        end do
      end do

    end do
  end do
  p2_fold = Half*p2_fold
  p2a_fold = Half*p2a_fold

end subroutine fold_two_pdm

subroutine ex1_init(k,n,ex1_table)

  use second_quantization, only: ex1, fase, lex_init, lex_next, lexrank

  integer(kind=iwp), intent(in) :: k, n
  type(ex1_struct), intent(out) :: ex1_table(:,:)
  integer(kind=iwp) :: counter, det, idet, my_ndet, p, pq, q, tmp

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

subroutine hdiag(h,g,diag,CSFDIA,nConf)

  use Constants, only: Zero

  integer(kind=iwp), intent(in):: nConf
  real(kind=wp) :: CSFDIA(nConf)
  real(kind=wp), intent(in) :: h(my_norb,my_norb), g(my_norb,my_norb,my_norb,my_norb)

  real(kind=wp), intent(out) :: diag(my_ndet)

  real(kind=wp), allocatable :: Da(:), Db(:)

  integer(kind=iwp), allocatable :: occa(:,:), occb(:,:)

  integer(kind=iwp) :: ia, ib, i, j, k, iaLow
! ToDo cleanup this!
! integer(kind=iwp) :: MinOpen, MaxOpen, nOpen, nAEl, nBEl, iDet, jDet, iCSF, iCNF, jCNF, iCSFOff
  integer(kind=iwp) :: MinOpen, MaxOpen, nOpen, nAEl,       iDet, jDet, iCSF, iCNF, jCNF, iCSFOff
  integer(kind=iwp) :: SumDet, SumCSF, ipat, ipat_match
  real(kind=wp) :: EAver, EDet

  integer(kind=iwp) :: occ(my_norb)

  call mma_allocate(Da,ndeta,label='Da')
  call mma_allocate(Db,ndetb,label='Db')

  call mma_allocate(occa,nela,ndeta,label='OccA')
  call mma_allocate(occb,nelb,ndetb,label='OccB')

  call build_occ_alpha(occa)
  call build_occ_beta (occb)
  !
  ! Alpha-string energies
  !
  do ia=1,ndeta
     Da(ia) = Zero
     do i=1,nela
        Da(ia) = Da(ia) + h(occa(i,ia),occa(i,ia))
     end do
     do i=1,nela
        do j=i+1,nela
           Da(ia) = Da(ia)                                     &
              + g(occa(i,ia),occa(i,ia),occa(j,ia),occa(j,ia)) &
              - g(occa(i,ia),occa(j,ia),occa(j,ia),occa(i,ia))
        end do
     end do
  end do
  !
  ! Beta-string energies
  !
  do ib=1,ndetb
     Db(ib) = Zero
     do i=1,nelb
        Db(ib) = Db(ib) + h(occb(i,ib),occb(i,ib))
     end do
     do i=1,nelb
        do j=i+1,nelb
           Db(ib) = Db(ib)                                     &
              + g(occb(i,ib),occb(i,ib),occb(j,ib),occb(j,ib)) &
              - g(occb(i,ib),occb(j,ib),occb(j,ib),occb(i,ib))
        end do
     end do
  end do
! Todo: Remove the computation of diag!(?)
  !
  ! Full determinant-product diagonal
  !
  k = 0
  do ib=1,ndetb
     iaLow = 1
     if (Mult==1) iaLow=ib
     do ia=iaLow,ndetA
        k = k + 1
!       Write (6,*) 'k=',k,' -> ib, ia:',ib,ia
        diag(k) = Da(ia) + Db(ib)
        do i=1,nela
           do j=1,nelb
              diag(k) = diag(k)                                   &
                 + g(occa(i,ia),occa(i,ia),                       &
                     occb(j,ib),occb(j,ib))
           end do
        end do
     end do
  end do

! Call RecPrt('Diag',' ',Diag,1,k)

! Construct the averaged determinant contribution of a configurtion and distribute it
! to all the CSFs of the same configuration.

  MinOpen=Mult-1
  MaxOpen=Min(my_nEl,my_norb)
  If (mod(MaxOpen-my_nEl,2) /= 0) MaxOpen=MaxOpen-1

! nOpen = Index of open-shell class. Equal the actual number of open shells.
  CSFDIA(:)=Zero
  iCSFOff=0
  sumDet = 0
  sumCSF = 0
  Do nOpen = MinOpen, MaxOpen, 2
     nAEl=(nOpen+(Mult-1))/2
!    nBEl=(nOpen-(Mult-1))/2
!    IDET = number of prototype determinants per configuration in this open-shell class.
!    ICSF = number of CSFs per configuration in this open-shell class.
!    ICNF = number of configurations in this open-shell class.

!
! Number of prototype determinants associated with this
! configuration (occupation pattern).  For singlets this is
! smaller than the full determinant count because spin-reversed
! determinant pairs are merged into a single prototype
! determinant.  Therefore use
!
!    iDet = nComb(ipat_match)
!
! rather than the full determinant multiplicity
!
!    binom(nOpen,nAEl)
!
! when constructing the configuration average.
!
     iDet = binom_coef(nAEl,nOpen)
     iCSF=ncsf_open_shells(nOpen,Mult-1)
     iCNF=nconf_open_shells(my_norb,my_nel,nOpen)
!    Write(6,'(A,5I6)') 'nOpen, nAEl, nBEl, nDet, nCSF=', nOpen, nAEl, nBEl, iDet, iCSF
     sumDet = sumDet + iDet*iCNF
     sumCSF = sumCSF + iCSF*iCNF
     Do jCNF=1, iCNF

        call configuration_from_index(nOpen,jCNF,occ)

        ipat_match = 0

        do ipat=1,npat
           if (all(occ(:) == occpat(:,ipat))) then
              ipat_match = ipat
              exit
           end if
        end do
!       write(u6,'(A,3I6,A,20I2)') 'nOpen,jCNF,ipat=',nOpen,jCNF,ipat_match,' occ=',occ(:)

        EAver = Zero

        iDet = nComb(ipat_match)
        do jDet=1,iDet
          call prototype_diagonal(ipat_match,jDet,eDet)
!         write(u6,*) jDet, eDet
          EAver = EAver + eDet
        end do
        EAver = EAver / real(iDet,wp)
!       write(u6,'(A,2I6,F20.10)') 'jCNF, ipat =', jCNF, ipat_match, EAver

!       write(u6,'(A,4I6)') 'nOpen,jCNF,ipat,nCSF=',  nOpen,jCNF,ipat_match,iCSF

        CSFDIA(iCSFOff+1:iCSFOff+iCSF) = EAver
        iCSFOff = iCSFOff + iCSF
     End Do
  End Do
! Write(6,'(A,2I6)') 'SumDet, SumCSF=', SumDet, SumCSF

! Call RecPrt('CSFDIA(Farolad)',' ',CSFDIA,1,nConf)

  call mma_deallocate(occb)
  call mma_deallocate(occa)
  call mma_deallocate(Db)
  call mma_deallocate(Da)

contains

subroutine prototype_diagonal(ipat,jDet,e)

   integer(kind=iwp), intent(in) :: ipat
   integer(kind=iwp), intent(in) :: jDet
   real(kind=wp), intent(out)    :: e

   integer(kind=iwp) :: nOpen
   integer(kind=iwp) :: nCombPat
   integer(kind=iwp), allocatable :: comb(:,:)

   nOpen = nsing(ipat)

   if (nOpen == 0) then

      call combination_diagonal(ipat,htu,gtuvx,e)

   else

      call mma_allocate(comb,nOpen,ncomb(ipat),label='Comb')

      call pattern_combinations(ipat,nOpen,nCombPat,comb)

      call combination_diagonal(ipat,htu,gtuvx,e,comb(:,jDet))

      call mma_deallocate(comb)

   end if

end subroutine


  subroutine build_occ_alpha(occ)

    integer(kind=iwp), intent(out) :: occ(nela,ndeta)

    integer(kind=iwp) :: ia, ipq, nocc

    do ia=1,ndeta

       nocc = 0

       do ipq=1,max_ex1a

          if (ex1_a(ipq,ia)%p == ex1_a(ipq,ia)%q .and. &
              ex1_a(ipq,ia)%rank == ia) then

             nocc = nocc + 1
             occ(nocc,ia) = ex1_a(ipq,ia)%p

          end if

       end do

       if (nocc /= nela) then
          write(u6,*) 'build_occ_alpha: wrong occupation count'
          call Abend()
       end if

    end do

  end subroutine build_occ_alpha

  subroutine build_occ_beta(occ)

    integer(kind=iwp), intent(out) :: occ(nelb,ndetb)

    integer(kind=iwp) :: ib, ipq, nocc

    do ib=1,ndetb

       nocc = 0

       do ipq=1,max_ex1b

          if (ex1_b(ipq,ib)%p == ex1_b(ipq,ib)%q .and. &
              ex1_b(ipq,ib)%rank == ib) then

             nocc = nocc + 1
             occ(nocc,ib) = ex1_b(ipq,ib)%p

          end if

       end do

       if (nocc /= nelb) then
          write(u6,*) 'build_occ_beta: wrong occupation count'
          call Abend()
       end if

    end do

  end subroutine build_occ_beta


end subroutine hdiag

subroutine build_pattern_diagonal_approx(h,g,diag)

! Occupation-pattern (Olsen-type) diagonal.
! Does NOT represent the exact CSF diagonal.

  use Constants, only: Zero

  real(kind=wp), intent(in) :: h(my_norb,my_norb), g(my_norb,my_norb,my_norb,my_norb)
  real(kind=wp), intent(out) :: diag(npat)

  integer(kind=iwp) :: ipat, p, q

  integer(kind=iwp) :: np, nq, exch

  do ipat=1,npat

     diag(ipat) = Zero

     !
     ! One-electron contribution
     !
     do p=1,my_norb

        np = occpat(p,ipat)

        diag(ipat) = diag(ipat) + real(np,wp)*h(p,p)

     end do

     !
     ! Two-electron contribution
     !
     do p=1,my_norb-1

        np = occpat(p,ipat)
        if (np == 0) cycle

        do q=p+1,my_norb

           nq = occpat(q,ipat)
           if (nq == 0) cycle

           !
           ! Coulomb term
           !
           diag(ipat) = diag(ipat) + real(np*nq,wp) * g(p,p,q,q)

           !
           ! Exchange coefficient
           !
           if (np == 2 .and. nq == 2) then
              exch = 2
           else
              exch = 1
           end if

           diag(ipat) = diag(ipat) - real(exch,wp) * g(p,q,q,p)

        end do

     end do

  end do

end subroutine build_pattern_diagonal_approx

subroutine det_occ_pattern(ia,ib,occ)

  integer(kind=iwp), intent(in) :: ia, ib

  integer(kind=iwp), intent(out) :: occ(my_norb)

  integer(kind=iwp) :: p, pq

  occ(:) = 0

  !
  ! Alpha occupations
  !

  do pq=1,max_ex1a

     if (ex1_a(pq,ia)%p == ex1_a(pq,ia)%q .and. &
         ex1_a(pq,ia)%rank == ia) then

        p = ex1_a(pq,ia)%p
        occ(p) = occ(p) + 1

     end if

  end do

  !
  ! Beta occupations
  !

  do pq=1,max_ex1b

     if (ex1_b(pq,ib)%p == ex1_b(pq,ib)%q .and. &
         ex1_b(pq,ib)%rank == ib) then

        p = ex1_b(pq,ib)%p
        occ(p) = occ(p) + 1

     end if

  end do

end subroutine det_occ_pattern

subroutine verify_occ_patterns()

  integer(kind=iwp) :: ia, ib, pq, p

  integer(kind=iwp) :: occ(my_norb)

  integer(kind=iwp) :: na, nb

  do ib=1,ndetb

     do ia=1,ndeta

        occ(:) = 0

        !
        ! Recover alpha occupations
        !

        na = 0

        do pq=1,max_ex1a

           if (ex1_a(pq,ia)%p == ex1_a(pq,ia)%q .and. &
               ex1_a(pq,ia)%rank == ia) then

              p = ex1_a(pq,ia)%p

              occ(p) = occ(p) + 1
              na = na + 1

           end if

        end do

        !
        ! Recover beta occupations
        !

        nb = 0

        do pq=1,max_ex1b

           if (ex1_b(pq,ib)%p == ex1_b(pq,ib)%q .and. &
               ex1_b(pq,ib)%rank == ib) then

              p = ex1_b(pq,ib)%p

              occ(p) = occ(p) + 1
              nb = nb + 1

           end if

        end do

        !
        ! Consistency checks
        !

        if (na /= nela) then
           write(u6,*) 'verify_occ_patterns: wrong alpha count'
           write(u6,*) 'ia=',ia,' na=',na,' nela=',nela
           call Abend()
        end if

        if (nb /= nelb) then
           write(u6,*) 'verify_occ_patterns: wrong beta count'
           write(u6,*) 'ib=',ib,' nb=',nb,' nelb=',nelb
           call Abend()
        end if

        if (sum(occ) /= my_nel) then
           write(u6,*) 'verify_occ_patterns: wrong electron count'
           write(u6,*) 'ia=',ia,' ib=',ib
           write(u6,*) 'sum(occ)=',sum(occ)
           write(u6,*) 'my_nel  =',my_nel
           call Abend()
        end if

        if (maxval(occ) > 2) then
           write(u6,*) 'verify_occ_patterns: occupation > 2'
           write(u6,*) 'ia=',ia,' ib=',ib
           call Abend()
        end if

     end do

  end do

end subroutine verify_occ_patterns

subroutine build_patterns()

  use stdalloc, only: mma_allocate, mma_deallocate

  integer(kind=iwp) :: ia, ib, k, ipat, jpat

  integer(kind=iwp) :: occ(my_norb)

  integer(kind=iwp), allocatable :: tmp(:,:)

  logical :: found

  !
  ! Worst case: every determinant product generates
  ! its own occupation pattern.
  !

  call mma_allocate(occpat,my_norb,my_ndet,label='OccPat')
  call mma_allocate(ipat_of_det,my_ndet,label='IPatDet')

  npat = 0
  k    = 0

  do ib=1,ndetb

     do ia=1,ndeta

        k = k + 1

        call det_occ_pattern(ia,ib,occ)

        found = .false.

        do jpat=1,npat

           if (all(occ(:) == occpat(:,jpat))) then

              ipat_of_det(k) = jpat
              found = .true.

              exit

           end if

        end do

        if (.not. found) then

           npat = npat + 1

           occpat(:,npat) = occ(:)

           ipat_of_det(k) = npat

        end if

     end do

  end do

  !
  ! Compress occpat to its actual size.
  !

  call mma_allocate(tmp,my_norb,npat,label='TmpPat')

  tmp(:,:) = occpat(:,1:npat)

  call mma_deallocate(occpat)

  call mma_allocate(occpat,my_norb,npat,label='OccPat')

  occpat(:,:) = tmp(:,:)

  call mma_deallocate(tmp)

  !
  ! Consistency checks
  !

  do k=1,my_ndet
     ipat = ipat_of_det(k)

     if (ipat < 1 .or. ipat > npat) then
        write(u6,*) 'build_patterns: invalid pattern index'
        write(u6,*) 'k=',k,' ipat=',ipat,' npat=',npat
        call Abend()
     end if

  end do

#ifdef _DEBUGPRINT_
  do k = 1, nPat
     write(u6,'(A,20I2)') 'occ=',occpat(:,k)
  end do
  write(u6,'(A,I10)') ' NDet = ', my_ndet
  write(u6,'(A,I10)') ' NPat = ', npat
  write(u6,'(A,F12.4)') ' Compression = ', real(my_ndet,wp)/real(npat,wp)
#endif

#ifdef _DEBUGPRINT_
do ipat=1,npat
   write(u6,'(A,I4,A,20I2)') &
        'ipat=',ipat, &
        ' occ=',occpat(:,ipat)
end do
#endif


end subroutine build_patterns

!#define _Print_Combinations_
subroutine analyse_patterns()

  integer(kind=iwp) :: ipat, p

#ifdef _Print_Combinations_
integer(kind=iwp) :: nopen, ncomb_pat, icomb
integer(kind=iwp), allocatable :: Comb(:,:)
integer(kind=iwp) :: occa_tmp(nela), occb_tmp(nelb)
#endif

  call mma_allocate(ndoub,npat,label='NDoub')
  call mma_allocate(nsing,npat,label='NSing')
  call mma_allocate(nspin_comb,npat,label='NSpinComb')
  call mma_allocate(ncomb,npat,label='NComb')

  do ipat=1,npat

     ndoub(ipat) = 0
     nsing(ipat) = 0

     do p=1,my_norb

        select case (occpat(p,ipat))

        case (2)

           ndoub(ipat) = ndoub(ipat) + 1

        case (1)

           nsing(ipat) = nsing(ipat) + 1

        end select

     end do

     if (2*ndoub(ipat)+nsing(ipat) /= my_nel) then

        write(u6,*) 'analyse_patterns: electron count error'
        write(u6,*) 'ipat=',ipat

        call Abend()

     end if
!
! Lucia uses spin combinations only for singlets (MS2 = 0).
! For all other spins the code falls back to the ordinary
! SD representation. Therefore nComb equals the number of
! alpha/beta assignments among the open shells.
!
     nspin_comb(ipat) = spin_comb_count(nsing(ipat),mult-1)
     if (mult /= 1) then
       !
       ! Lucia uses spin combinations only for singlets.
       ! For non-singlets the code falls back to the ordinary
       ! SD representation
       !
       ncomb(ipat) = nspin_comb(ipat)
     else
       if (nsing(ipat) == 0) then
         ncomb(ipat) = 1
       else
         ncomb(ipat) = nspin_comb(ipat)/2
       end if
     end if

#ifdef _Print_Combinations_
   call mma_allocate(comb,nsing(ipat),nspin_comb(ipat),label='Comb')
   call pattern_combinations(ipat,nopen,ncomb_pat,comb)

   do icomb=1,ncomb_pat

      call combination_occupations(ipat,comb(:,icomb),occa_tmp,occb_tmp)

      write(u6,'(A,I3)') 'Combination ',icomb
      write(u6,'(A,20I3)') 'Alpha:',occa_tmp
      write(u6,'(A,20I3)') 'Beta :',occb_tmp

   end do

   call mma_deallocate(comb)
#endif

  end do

#ifdef _DEGUGPRINT_
write(u6,*)
write(u6,*) 'IPAT  NOPEN  NDETALL  NDETPROTO  NCSF'

do ipat=1,npat

   write(u6,'(5I8)') &
        ipat, &
        nsing(ipat), &
        count(ipat_of_det == ipat), &
        ncomb(ipat), &
        ncsf_open_shells(nsing(ipat),mult-1)

end do
#endif

call mma_allocate(ibcomb,npat,label='IBComb')

ibcomb(1) = 1

do ipat=2,npat
   ibcomb(ipat) = ibcomb(ipat-1) + ncomb(ipat-1)
end do

if (mult == 1) then
   ncomb_tot = ndeta*(ndeta+1)/2
else
   ncomb_tot = ndeta*ndetb
end if

call build_ictsdt()


end subroutine analyse_patterns

integer(kind=iwp) function spin_comb_count(iopen,ms2)

  integer(kind=iwp), intent(in) :: iopen, ms2
  integer(kind=iwp) :: iael, ibel

  iael = (iopen+ms2)/2
  ibel = (iopen-ms2)/2

  if (iael < 0 .or. ibel < 0) then
     spin_comb_count = 0
  else if (iael+ibel /= iopen) then
     spin_comb_count = 0
  else
     spin_comb_count = binom_coef(iael,iopen)
  end if

end function spin_comb_count

subroutine spncom_faroald(nopen,ms2,ncomb,comb)

  integer(kind=iwp), intent(in) :: nopen, ms2
  integer(kind=iwp), intent(in) :: ncomb
  integer(kind=iwp), intent(out) ::  comb(nopen,ncomb)
  integer(kind=iwp) :: i,j, add, nalpha, icomb
  integer(kind=iwp) :: work(nopen)
  integer(kind=iwp) :: mx

  work(:) = 0
  mx = 2**nopen
  icomb = 0

  do i=1,mx
     if (i > 1) then
        add = 1
        j   = 0
        do while (add == 1)
           j = j + 1
           if (work(j) == 1) then
              work(j) = 0
           else
              work(j) = 1
              add = 0
           end if
        end do
     end if

     nalpha = sum(work)

     if (2*nalpha-nopen == ms2) then
       if (mult /= 1 .or. work(1) == 1) then
         icomb = icomb + 1
         comb(:,icomb) = work(:)
       end if
     end if

  end do

end subroutine spncom_faroald


subroutine pattern_combinations(ipat,nopen,ncomb_pat,comb)

  integer(kind=iwp), intent(in) :: ipat
  integer(kind=iwp), intent(out) :: nopen, ncomb_pat
  integer(kind=iwp), intent(out) :: comb(:,:)

  integer(kind=iwp) :: p

  nopen = 0

  do p=1,my_norb
     if (occpat(p,ipat) == 1) nopen = nopen + 1
  end do

  ncomb_pat = ncomb(ipat)

  call spncom_faroald(nopen,mult-1,ncomb_pat,comb)

end subroutine pattern_combinations

subroutine combination_occupations(ipat,comb,occa,occb)

  integer(kind=iwp), intent(in) :: ipat
  integer(kind=iwp), intent(in), optional :: comb(:)

  integer(kind=iwp), intent(out) :: occa(nela)
  integer(kind=iwp), intent(out) :: occb(nelb)

  integer(kind=iwp) :: ia
  integer(kind=iwp) :: ib
  integer(kind=iwp) :: ic
  integer(kind=iwp) :: p

  ia = 0
  ib = 0
  ic = 0

  do p=1,my_norb

     select case (occpat(p,ipat))

     case (2)

        ia = ia + 1
        ib = ib + 1

        occa(ia) = p
        occb(ib) = p

     case (1)

        ic = ic + 1

        if (.not. present(comb)) then
          write (u6,*) 'combination_occupations. comb missing'
          call Abend()
        else if (comb(ic) == 1) then

           ia = ia + 1
           occa(ia) = p

        else

           ib = ib + 1
           occb(ib) = p

        end if

     end select

  end do

end subroutine combination_occupations

subroutine combination_diagonal(ipat,h,g,e,comb)

  real(kind=wp), intent(in) :: h(my_norb,my_norb)
  real(kind=wp), intent(in) :: g(my_norb,my_norb,my_norb,my_norb)

  integer(kind=iwp), intent(in) :: ipat
  integer(kind=iwp), intent(in), optional :: comb(:)

  real(kind=wp), intent(out) :: e

  integer(kind=iwp) :: occa(nela)
  integer(kind=iwp) :: occb(nelb)

  integer(kind=iwp) :: i,j

  call combination_occupations(ipat,comb,occa,occb)

  e = Zero

  do i=1,nela
     e = e + h(occa(i),occa(i))
  end do

  do i=1,nelb
     e = e + h(occb(i),occb(i))
  end do

  do i=1,nela
     do j=i+1,nela
        e = e + g(occa(i),occa(i),occa(j),occa(j))
        e = e - g(occa(i),occa(j),occa(j),occa(i))
     end do
  end do

  do i=1,nelb
     do j=i+1,nelb
        e = e + g(occb(i),occb(i),occb(j),occb(j))
        e = e - g(occb(i),occb(j),occb(j),occb(i))
     end do
  end do

  do i=1,nela
     do j=1,nelb
        e = e + g(occa(i),occa(i),occb(j),occb(j))
     end do
  end do

end subroutine combination_diagonal

subroutine build_pattern_diagonal(h,g,diag)

  use Constants, only: Zero

  real(kind=wp), intent(in) :: h(my_norb,my_norb), g(my_norb,my_norb,my_norb,my_norb)
  real(kind=wp), intent(out) :: diag(npat)

  integer(kind=iwp) :: ipat, p,q, np,nq, exch


  do ipat=1,npat
     diag(ipat) = Zero

     !
     ! One-electron contribution
     !

     do p=1,my_norb
        np = occpat(p,ipat)
        diag(ipat) = diag(ipat) + real(np,wp)*h(p,p)
     end do

     !
     ! Two-electron contribution
     !

     do p=1,my_norb-1
        np = occpat(p,ipat)
        if (np == 0) cycle
        do q=p+1,my_norb
           nq = occpat(q,ipat)
           if (nq == 0) cycle
           !
           ! Coulomb
           !
           diag(ipat) = diag(ipat) + real(np*nq,wp) * g(p,p,q,q)
           !
           ! Exchange
           !
           if (np == 2 .and. nq == 2) then
              exch = 2
           else
              exch = 1
           end if
           diag(ipat) = diag(ipat) - real(exch,wp) * g(p,q,q,p)

        end do
     end do
     !
     ! On-site double occupations
     !

     do p=1,my_norb
        if (occpat(p,ipat) == 2) then
           diag(ipat) = diag(ipat) + g(p,p,p,p)
        end if
     end do
  end do

end subroutine build_pattern_diagonal

subroutine build_ictsdt()

  integer(kind=iwp) :: ia
  integer(kind=iwp) :: ib
  integer(kind=iwp) :: iaLow
  integer(kind=iwp) :: k

  call mma_allocate(ictsdt,ncomb_tot,label='ICTSDT')

  k = 0

  do ib=1,ndetb

     iaLow = 1
     if (mult == 1) iaLow = ib

     do ia=iaLow,ndeta

        k = k + 1

        ictsdt(k) = k

     end do

  end do

end subroutine build_ictsdt


integer(kind=iwp) function nconf_open_shells(norb,nel,nopen)

  integer(kind=iwp), intent(in) :: norb
  integer(kind=iwp), intent(in) :: nel
  integer(kind=iwp), intent(in) :: nopen

  integer(kind=iwp) :: ndoub

  if (mod(nel-nopen,2) /= 0) then
     nconf_open_shells = 0
     return
  end if

  ndoub = (nel-nopen)/2

  nconf_open_shells = &
       binom_coef(nopen,norb) * &
       binom_coef(ndoub,norb-nopen)

end function nconf_open_shells

integer(kind=iwp) function ncsf_open_shells(nopen,ms2)

  integer(kind=iwp), intent(in) :: nopen
  integer(kind=iwp), intent(in) :: ms2

  integer(kind=iwp) :: k

  ! parity check
  if (mod(nopen-ms2,2) /= 0) then
     ncsf_open_shells = 0
     return
  end if

  k = (nopen-ms2)/2

  ncsf_open_shells = binom_coef(k,nopen)

  if (k > 0) then
     ncsf_open_shells = ncsf_open_shells - &
                        binom_coef(k-1,nopen)
  end if

end function ncsf_open_shells


subroutine configuration_from_index(nOpen,jCNF,occ)

  integer(kind=iwp), intent(in)  :: nOpen
  integer(kind=iwp), intent(in)  :: jCNF
  integer(kind=iwp), intent(out) :: occ(my_norb)

  integer(kind=iwp) :: iconf_exp(my_nel)
  integer(kind=iwp) :: occ_local(my_norb)

  integer(kind=iwp) :: i
  integer(kind=iwp) :: nopen_local
  integer(kind=iwp) :: count_conf
  integer(kind=iwp) :: nonew

  occ(:) = 0

  !
  ! Initial configuration
  !
  call nxt_conf_faroald(iconf_exp,my_nel,my_norb,1,nonew)

  count_conf = 0

  do while (nonew == 0)

     !
     ! Expanded occupation list -> occupation pattern
     !
     occ_local(:) = 0

     do i=1,my_nel
        occ_local(iconf_exp(i)) = occ_local(iconf_exp(i)) + 1
     end do

     !
     ! Lucia only allows occupations 0,1,2
     !
     if (maxval(occ_local) <= 2) then

        nopen_local = count(occ_local(:) == 1)

        if (nopen_local == nOpen) then

           count_conf = count_conf + 1

           if (count_conf == jCNF) then

              occ(:) = occ_local(:)

              return

           end if

        end if

     end if

     !
     ! Next configuration in Lucia order
     !
     call nxt_conf_faroald(iconf_exp,my_nel,my_norb,0,nonew)

  end do

  write(u6,*) 'configuration_from_index: jCNF out of range'
  write(u6,*) 'nOpen=',nOpen,' jCNF=',jCNF
  call Abend()

end subroutine configuration_from_index

subroutine nxt_conf_faroald(iconf, nel, norb, ini, nonew)

  integer(kind=iwp), intent(in)    :: nel
  integer(kind=iwp), intent(in)    :: norb
  integer(kind=iwp), intent(in)    :: ini

  integer(kind=iwp), intent(inout) :: iconf(nel)

  integer(kind=iwp), intent(out)   :: nonew

  integer(kind=iwp) :: i
  integer(kind=iwp) :: iadd
  integer(kind=iwp) :: iel
  integer(kind=iwp) :: increase
  integer(kind=iwp) :: jorb
  integer(kind=iwp) :: ndouble

  if (ini == 1) then

     if (nel <= 2*norb) then

        nonew = 0

        ndouble = nel/2

        do i=1,ndouble
           iconf(2*i-1) = i
           iconf(2*i  ) = i
        end do

        if (2*ndouble /= nel) then
           iconf(2*ndouble+1) = ndouble+1
        end if

     else

        nonew = 1

     end if

  else

     iadd = 1
     iel  = 0

     do while (iadd == 1)

        iel = iel + 1

        increase = 0

        if (iel < nel) then

           if (iconf(iel) < iconf(iel+1)-1) then

              increase = 1

           else if (iconf(iel) == iconf(iel+1)-1) then

              if (iel == nel-1) then

                 increase = 1

              else if (iconf(iel+1) /= iconf(iel+2)) then

                 increase = 1

              end if

           end if

        else if (iel == nel) then

           if (iconf(iel) < norb) then

              increase = 1

           else

              nonew = 1
              exit

           end if

        else

           nonew = 1
           exit

        end if

        if (increase == 1) then

           nonew = 0

           iconf(iel) = iconf(iel) + 1

           ndouble = (iel-1)/2

           do jorb=1,ndouble
              iconf(2*jorb-1) = jorb
              iconf(2*jorb  ) = jorb
           end do

           if (2*ndouble < iel-1) then
              iconf(iel-1) = ndouble + 1
           end if

           iadd = 0

        end if

     end do

  end if

end subroutine nxt_conf_faroald

subroutine LRs_init(p,q,my_nel,my_norb,L,R,sgn,counter)
! for a pair of orbitals p and q, and determinants
! generated by my_nel electrons in my_norb spin orbitals,
! re-enumerate all non-vanishing couples connected
! through: jdet = E_pq idet. The number of couples
! is given by n_det, L(i) = jdet and R(i) = idet,
! where i is a counter from 1 to n_det. The sign of
! jdet is stored in the sgn array.

  use second_quantization, only: ex1, fase, lex_init, lex_next, lexrank

  integer(kind=iwp), intent(in) :: p, q, my_nel, my_norb
  integer(kind=iwp), intent(out) :: L(:), R(:), sgn(:), counter
  integer(kind=iwp) :: det, idet, my_ndet, tmp

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

subroutine Faroald_epq_psi(p,q,psi,sgm)

  integer(kind=iwp), intent(in) :: p,q
  real(kind=wp), intent(in)  :: psi(ndeta,ndetb)
  real(kind=wp), intent(out) :: sgm(ndeta,ndetb)

  integer(kind=iwp) :: ja,jb
  integer(kind=iwp) :: ia,ib
  integer(kind=iwp) :: k

  sgm = Zero

  !
  ! alpha excitation part
  !
  do ja=1,ndeta
    do k=1,max_ex1a

      if (ex1_a(k,ja)%p /= p) cycle
      if (ex1_a(k,ja)%q /= q) cycle

      ia = ex1_a(k,ja)%rank

      sgm(ia,:) = sgm(ia,:) +         &
           real(ex1_a(k,ja)%sgn,wp) * &
           psi(ja,:)

    end do
  end do

  !
  ! beta excitation part
  !
  do jb=1,ndetb
    do k=1,max_ex1b

      if (ex1_b(k,jb)%p /= p) cycle
      if (ex1_b(k,jb)%q /= q) cycle

      ib = ex1_b(k,jb)%rank

      sgm(:,ib) = sgm(:,ib) +         &
           real(ex1_b(k,jb)%sgn,wp) * &
           psi(:,jb)

    end do
  end do

end subroutine Faroald_epq_psi

subroutine sigma_update_1e(h1,psi,sigma)

  real(kind=wp), intent(in)  :: h1(:,:)
  real(kind=wp), intent(in)  :: psi(ndeta,ndetb)
  real(kind=wp), intent(out) :: sigma(ndeta,ndetb)

  integer(kind=iwp) :: p,q
  real(kind=wp) :: hpq
  real(kind=wp), allocatable :: work(:,:)
  real(kind=wp), parameter :: thr=1.0e-12_wp

  Call mma_allocate(Work,nDetA,nDetB,Label='Work')
  sigma(:,:) = Zero

  do p=1,my_norb
    do q=1,my_norb

      hpq = h1(p,q)

      if (abs(hpq) < thr) cycle

      call Faroald_Epq_Psi(p,q,psi,work)

      sigma(:,:) = sigma(:,:) + hpq*work(:,:)

    end do
  end do
  Call mma_deallocate(Work)

end subroutine sigma_update_1e

subroutine traci_faroald(T,nT,psi)

   integer(iwp), intent(in) :: nT
   real(wp), intent(in)    :: T(nT)
   real(wp), intent(inout) :: psi(ndeta,ndetb)

   real(wp), allocatable :: h1(:,:)
   real(wp), allocatable :: psi0(:,:), psi1(:,:), psi2(:,:)

   integer(iwp) :: k
   real(wp) :: tkk

   call mma_allocate(h1,my_norb,my_norb,Label='h1')
   call mma_allocate(psi0,ndeta,ndetb,Label='Psi0')
   call mma_allocate(psi1,ndeta,ndetb,Label='Psi1')
   call mma_allocate(psi2,ndeta,ndetb,Label='Psi2')

   Call RecPrt('Psi(Faroald)',' ',Psi,nDeta,nDetb)
   Call RecPrt('T',' ',T,my_nOrb,my_nOrb)
   do k=1,my_norb
      Write (6,*) 'k=',k
      call build_trow_operator(T,k,h1,tkk)
      Write (6,*) 'tkk=',tkk
      Call RecPrt('h1',' ',h1,my_nOrb,my_nOrb)

      psi0(:,:) = psi(:,:)

      call multiply_by_nk(tkk,k,psi0)
      Write (6,*) 'After multiply_by_nk'
      Call RecPrt('Psi0',' ',Psi0,nDetA,nDetb)

      call sigma_update_1e(h1,psi0,psi1)
      Write (6,*) 'After first sigma_update_1e'
      Call RecPrt('Psi1',' ',Psi1,nDetA,nDetb)

      call sigma_update_1e(h1,psi1,psi2)
      Write (6,*) 'After second sigma_update_1e'
      Call RecPrt('Psi2',' ',Psi2,nDetA,nDetb)

      psi (:,:)= psi0(:,:)+ psi1(:,:)+ Half*psi2(:,:)
      Call RecPrt('psi0(:,:)+ psi1(:,:)+ Half*psi2(:,:)',' ',Psi,nDetA,nDetb)

   end do

   call mma_deallocate(psi2)
   call mma_deallocate(psi1)
   call mma_deallocate(psi0)
   call mma_deallocate(h1)

end subroutine traci_faroald

subroutine build_trow_operator(T,k,h1,tkk)

   real(wp), intent(in)  :: T(:)
   integer(iwp), intent(in) :: k
   real(wp), intent(out) :: h1(my_norb,my_norb)
   real(wp), intent(out) :: tkk

   integer(iwp) :: p

   h1(:,:) = Zero

   tkk = T((k-1)*my_norb + k)

   do p=1,my_norb
      if (p == k) cycle
      h1(p,k) = T((k-1)*my_norb + p) / tkk
   end do

end subroutine

subroutine multiply_by_nk(tkk,k,psi)

   use Constants,   only: One
   use Definitions, only: wp, iwp

   real(kind=wp), intent(in)    :: tkk
   integer(kind=iwp), intent(in) :: k
   real(kind=wp), intent(inout) :: psi(ndeta,ndetb)

   integer(kind=iwp) :: ia, ib, pq
   integer(kind=iwp), allocatable :: occa(:), occb(:)

   real(kind=wp) :: fac(0:2)

   fac(0)=One
   fac(1)=tkk
   fac(2)=tkk*tkk

   call mma_allocate(occa,ndeta,label='OccA')
   call mma_allocate(occb,ndetb,label='OccB')

   occa(:)=0
   occb(:)=0

   do ia=1,ndeta
      do pq=1,max_ex1a
         if (ex1_a(pq,ia)%p == k .and. &
             ex1_a(pq,ia)%q == k .and. &
             ex1_a(pq,ia)%rank == ia) then
            occa(ia)=1
            exit
         end if
      end do
   end do

   do ib=1,ndetb
      do pq=1,max_ex1b
         if (ex1_b(pq,ib)%p == k .and. &
             ex1_b(pq,ib)%q == k .and. &
             ex1_b(pq,ib)%rank == ib) then
            occb(ib)=1
            exit
         end if
      end do
   end do

   do ib=1,ndetb
      do ia=1,ndeta
         psi(ia,ib)=fac(occa(ia)+occb(ib))*psi(ia,ib)
      end do
   end do

   call mma_deallocate(occb)
   call mma_deallocate(occa)

end subroutine multiply_by_nk

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
