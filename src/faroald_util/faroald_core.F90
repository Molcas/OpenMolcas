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

module Faroald_Core
! written by Steven Vancoillie, summer 2014
!
! The Faroald module handles sigma updates as the result of acting
! with the hamiltonian operator on a CI vector: s = H c, where s and c
! are CI expansions in determinant basis.
!
! The implementation follows the minimum operation count algorithm
! published by Olsen & Co in J. Chem. Phys. 89, 2185 (1988).

use Faroald_SecQuant, only: Faro_Binom
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, Half
use Definitions, only: wp, iwp, u6
#ifdef _PROF_
use, intrinsic :: iso_fortran_env, only: int64
use Definitions, only: u6
#endif

implicit none
private

! wavefunction info (set by Faro_Init; read-only outside this module)
integer(kind=iwp), protected :: Faro_nEl = 0, Faro_nOrb = 0, Faro_mult = 0, Faro_nela = 0, Faro_nelb = 0, Faro_nhoa = 0, &
                                Faro_nhob = 0, Faro_ndeta = 0, Faro_ndetb = 0, Faro_nDet = 0
! the tables of the full string/determinant spaces are set up (Faro_Init)
logical(kind=iwp), protected :: Faro_Full = .false.

! excitation tables
type ex1_struct
  integer(kind=iwp) :: p, q, sgn, rank
end type

type(ex1_struct), allocatable :: Faro_ex1_a(:,:), Faro_ex1_b(:,:)
integer(kind=iwp) :: Faro_max_ex1a, Faro_max_ex1b, Faro_max_ex2a, Faro_max_ex2b, Faro_max_LRs

#ifdef _PROF_
integer(kind=int64) :: nflop
#endif

integer(kind=iwp) :: Faro_npat

integer(kind=iwp), allocatable :: Faro_occpat(:,:),  &    ! (Faro_nOrb,Faro_npat)
                                  Faro_ipat_of_det(:)
integer(kind=iwp), allocatable :: Faro_ndoub(:), Faro_nsing(:)
!
! Faro_nSpin_Comb : raw alpha/beta assignments among open shells
! Faro_nComb      : Lucia-compatible counting of diagonal states
!
integer(kind=iwp), allocatable :: Faro_nSpin_Comb(:), Faro_nComb(:)
integer(kind=iwp), allocatable :: Faro_ibComb(:)
integer(kind=iwp) :: Faro_nComb_tot
integer(kind=iwp), allocatable :: Faro_ictsdt(:)

public :: Faro_ex1_a, Faro_ex1_b, Faro_ex1_init, Faro_Full, Faro_Free, Faro_Init, Faro_max_ex1a, &
  Faro_max_ex1b, Faro_max_ex2a, Faro_max_ex2b, Faro_max_LRs, Faro_mult, &
          Faro_nDet, Faro_nEl, Faro_nOrb, Faro_ndeta, Faro_ndetb, Faro_nela, Faro_nelb, Faro_nhoa, Faro_nhob, &
            Faro_sigma_update, Faro_hDiag, Faro_verify_occ_patterns, Faro_build_patterns, &
            Faro_analyse_patterns, Faro_ndoub, Faro_nsing, &
          Faro_build_pattern_diagonal, Faro_build_pattern_diagonal_approx, Faro_nSpin_Comb, Faro_nComb, Faro_ibComb, &
            Faro_nComb_tot, Faro_pattern_combinations
public :: Faro_npat, Faro_occpat, Faro_ipat_of_det
public :: Faro_combination_occupations
public :: Faro_combination_diagonal
public :: Faro_build_ictsdt, Faro_ictsdt

! Extensions to mma interfaces

interface mma_allocate
  module procedure :: ex1_mma_allo_2D, ex1_mma_allo_2D_lim
end interface
interface mma_deallocate
  module procedure :: ex1_mma_free_2D
end interface

public :: mma_allocate, mma_deallocate

contains

subroutine Faro_Init(nactel,nasht,ispin,full)
! Fills in the global data of the module, in preparation for the
! Faro_sigma_update routine.
! full: also the excitation tables and occupation patterns of the full
! (CAS) string and determinant spaces (Faro_sigma_update, densities, traci);
! not needed for the string-driven RAS code (Faroald_RAS), which only uses
! the sizes, the rank tables and the integrals.

use Faroald_SecQuant, only: Faro_rank_init

integer(kind=iwp), intent(in) :: nactel, nasht, ispin
logical(kind=iwp), intent(in) :: full

Faro_nEl = nactel
Faro_nOrb = nasht
Faro_mult = ispin

! compute alpha/beta subsets
Faro_nela = (Faro_nEl+(Faro_mult-1))/2
Faro_nhoa = Faro_nOrb-Faro_nela
Faro_nelb = (Faro_nEl-(Faro_mult-1))/2
Faro_nhob = Faro_nOrb-Faro_nelb
Faro_ndeta = Faro_Binom(Faro_nela,Faro_nOrb)
Faro_ndetb = Faro_Binom(Faro_nelb,Faro_nOrb)
Faro_nDet = Faro_ndeta*Faro_ndetb

! initialize binomial tables needed for ranking
call Faro_rank_init()

! number of possible non-destructive single excitations
Faro_max_ex1a = Faro_nela*(Faro_nhoa+1)
Faro_max_ex1b = Faro_nelb*(Faro_nhob+1)

! The number of non-vanishing elements in F is the number of
! non-vanashing single and double excitations from a string, and
! computed as the number of strings which differ in two, one, or zero
! spin orbitals from a given string.
Faro_max_ex2a = (Faro_nela*(Faro_nela-1)*(Faro_nhoa)*(Faro_nhoa-1))/4+Faro_nela*Faro_nhoa+1
Faro_max_ex2b = (Faro_nelb*(Faro_nelb-1)*(Faro_nhob)*(Faro_nhob-1))/4+Faro_nelb*Faro_nhob+1

Faro_Full = full
if (full) then

! For SG_Epq_Psi, we use an excitation table for the beta strings.
call mma_allocate(Faro_ex1_b,Faro_max_ex1b,Faro_ndetb,label='ex1_b')
call Faro_ex1_init(Faro_nelb,Faro_nOrb,Faro_ex1_b)

! For sigma2, we use an excitation table for the alpha strings, unless
! Ms = 0 (Faro_mult = 1), as for singlet spins sigma2 is not computed and
! just taken as the transpose of SG_Epq_Psi.
!if (Faro_mult /= 1) then
call mma_allocate(Faro_ex1_a,Faro_max_ex1a,Faro_ndeta,label='ex1_a')
call Faro_ex1_init(Faro_nela,Faro_nOrb,Faro_ex1_a)
!end if

call Faro_verify_occ_patterns()
call Faro_build_patterns()
call Faro_analyse_patterns()

end if

! For sigma3, we need to construct a list of determinant couples L(i)
! = E_pq R(i) for a given p and q for the alpha strings. This is done
! by a subroutine that takes L, R, and sgn arrays. The maximum size of
! these arrays is when p=q, at which point the number of couples is:
Faro_max_LRs = Faro_Binom(Faro_nela-1,Faro_nOrb-1)


end subroutine Faro_Init


subroutine Faro_Free()
! Frees the global data of the module.


call mma_deallocate(Faro_ex1_a,safe='*')
call mma_deallocate(Faro_ex1_b,safe='*')
call mma_deallocate(Faro_ipat_of_det,safe='*')
call mma_deallocate(Faro_occpat,safe='*')
call mma_deallocate(Faro_ndoub,safe='*')
call mma_deallocate(Faro_nsing,safe='*')
call mma_deallocate(Faro_nSpin_Comb,safe='*')
call mma_deallocate(Faro_nComb,safe='*')
call mma_deallocate(Faro_ibComb,safe='*')
call mma_deallocate(Faro_ictsdt,safe='*')

end subroutine Faro_Free



subroutine Faro_sigma_update(h,g,sgm,psi)
! The sigma update routine performs the following operation:
! |sgm> = H |psi>, with the hamiltonian defined by
! H = sum_tu h(t,u) E_tu + sum_tuvx g(t,u,v,x) E_tuvx.

  ! integrals
  real(kind=wp), intent(in) :: h(Faro_nOrb,Faro_nOrb), g(Faro_nOrb,Faro_nOrb,Faro_nOrb,Faro_nOrb)
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
  call par_range(Faro_ndeta,iasta,iaend)
  call par_range(Faro_ndetb,ibsta,ibend)

# ifdef _PROF_
  ! initialize flop count
  nflop = 0

  call timing(t1_cpu,tot_cpu,t1_wall,tot_wall)
# endif

  ! Set the sigma vector to 0
  sgm = Zero

  ! First, construct a new effective one-electron integral matrix:
  ! k_tu = h_tu - 1/2 sum_v g_tvvu, to be used with sigma1/sigma2.
  call mma_allocate(k,Faro_nOrb,Faro_nOrb,label='k')
  do u=1,Faro_nOrb
    do t=1,Faro_nOrb
      ! g_tvvu = g_vtvu
      k(t,u) = Zero
      do v=1,Faro_nOrb
        k(t,u) = k(t,u)+g(v,t,v,u)
      end do
      k(t,u) = h(t,u)-Half*k(t,u)
    end do
  end do

  ! Second, for sigma2 and sigma3, we are better off with the transpose
  ! of sgm and/or psi, so allocate and assign them here.
  call mma_allocate(psiT,Faro_ndetb,Faro_ndeta,label='psiT')
  call dtrans(Faro_ndeta,Faro_ndetb,psi,Faro_ndeta,psiT,Faro_ndetb)

  ! Now the actual contributions to sigma are computed. For a singlet
  ! (Faro_mult = 1), sigma2 is not computed and sigma3 will only do half the
  ! work.  But at the end, the transpose of sigma has to be added to
  ! sigma. This should be more efficient than computing everything, but
  ! note that the transpose operation is also included in the timings
  ! used to compute the flop efficiency.

  call sigma1(k,g,sgm,psi,ibsta,ibend)

  if (Faro_mult /= 1) then
    ! we need efficient access to sgm by using the transpose
    call mma_allocate(sgmT,Faro_ndetb,Faro_ndeta,label='sgmT')
    call dtrans(Faro_ndeta,Faro_ndetb,sgm,Faro_ndeta,sgmT,Faro_ndetb)

    call sigma2(k,g,sgmT,psiT,iasta,iaend)

    call dtrans(Faro_ndetb,Faro_ndeta,sgmT,Faro_ndetb,sgm,Faro_ndeta)
    call mma_deallocate(sgmT)
  end if

  call sigma3(g,sgm,psiT,ibsta,ibend)

  ! sum over all processes
  call gadgop(sgm,Faro_ndeta*Faro_ndetb,'+')

  if (Faro_mult == 1) then
    ! for Ms = 0 (only used for singlet), sgm := sgm + sgm^T
    call transadd(Faro_ndeta,sgm,Faro_ndeta)
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

end subroutine Faro_sigma_update

subroutine sigma1(k,g,sgm,psi,ibsta,ibend)
! sigma1 = sum_jb sum_tu <jb|E_tu|ib> (h_kl - 1/2 sum_v <tv|vx>) C(ia,jb)
!    + 1/2 sum_jb sum_tuvx <jb|E_tu E_vx|ib> g_tuvx C(ia,jb)

  ! integrals
  real(kind=wp), intent(in) :: k(Faro_nOrb,Faro_nOrb), g(Faro_nOrb,Faro_nOrb,Faro_nOrb,Faro_nOrb)
  ! wavefunctions
  real(kind=wp), intent(inout) :: sgm(:,:)
  real(kind=wp), intent(in) :: psi(:,:)
  integer(kind=iwp), intent(in) :: ibsta, ibend
  ! local variables
  integer(kind=iwp) :: ib, jb, kb, t, u, v, x, tu, vx, sgn_tu, sgn_vx
  real(kind=wp), allocatable :: f(:)

  call mma_allocate(f,Faro_ndetb,label='f')

  do ib=ibsta,ibend
    ! f array construction
    f = Zero
    do tu=1,Faro_max_ex1b
      t = Faro_ex1_b(tu,ib)%p
      u = Faro_ex1_b(tu,ib)%q
      sgn_tu = Faro_ex1_b(tu,ib)%sgn
      kb = Faro_ex1_b(tu,ib)%rank
      f(kb) = f(kb)+sgn_tu*k(t,u)
      do vx=1,Faro_max_ex1b
        v = Faro_ex1_b(vx,kb)%p
        x = Faro_ex1_b(vx,kb)%q
        sgn_vx = Faro_ex1_b(vx,kb)%sgn
        jb = Faro_ex1_b(vx,kb)%rank
        f(jb) = f(jb)+Half*sgn_tu*sgn_vx*g(v,x,t,u)
      end do
    end do
    ! sigma addition
    kb = 0
    do jb=1,Faro_ndetb
      if (f(jb) /= Zero) then
        kb = kb+1
#       ifdef _PROF_
        nflop = nflop+2*Faro_ndeta
#       endif
        sgm(1:Faro_ndeta,ib) = sgm(1:Faro_ndeta,ib)+f(jb)*psi(1:Faro_ndeta,jb)
      end if
    end do
    if (kb > Faro_max_ex2b) stop 'exceeded max double excitations'
  end do

  call mma_deallocate(f)

end subroutine sigma1

subroutine sigma2(k,g,sgm,psi,iasta,iaend)
! sigma2 = sum_ja sum_tu <ja|E_tu|ia> (h_kl - 1/2 sum_v <tv|vx>) C(ja,ib)
!    + 1/2 sum_ja sum_tuvx <ja|E_tu E_vx|ia> g_tuvx C(ja,ib)

  ! integrals
  real(kind=wp), intent(in) :: k(Faro_nOrb,Faro_nOrb), g(Faro_nOrb,Faro_nOrb,Faro_nOrb,Faro_nOrb)
  ! wavefunctions
  real(kind=wp), intent(inout) :: sgm(:,:)
  real(kind=wp), intent(in) :: psi(:,:)
  integer(kind=iwp), intent(in) :: iasta, iaend
  ! local variables
  integer(kind=iwp) :: ia, ja, ka, t, u, v, x, tu, vx, sgn_tu, sgn_vx
  real(kind=wp), allocatable :: f(:)

  call mma_allocate(f,Faro_ndeta,label='f')

  do ia=iasta,iaend
    ! f array construction
    f = Zero
    do tu=1,Faro_max_ex1a
      t = Faro_ex1_a(tu,ia)%p
      u = Faro_ex1_a(tu,ia)%q
      sgn_tu = Faro_ex1_a(tu,ia)%sgn
      ka = Faro_ex1_a(tu,ia)%rank
      f(ka) = f(ka)+sgn_tu*k(t,u)
      do vx=1,Faro_max_ex1a
        v = Faro_ex1_a(vx,ka)%p
        x = Faro_ex1_a(vx,ka)%q
        sgn_vx = Faro_ex1_a(vx,ka)%sgn
        ja = Faro_ex1_a(vx,ka)%rank
        f(ja) = f(ja)+Half*sgn_tu*sgn_vx*g(v,x,t,u)
      end do
    end do
    ! sigma addition
    ka = 0
    do ja=1,Faro_ndeta
      if (f(ja) /= Zero) then
        ka = ka+1
#       ifdef _PROF_
        nflop = nflop+2*Faro_ndeta
#       endif
        sgm(1:Faro_ndetb,ia) = sgm(1:Faro_ndetb,ia)+f(ja)*psi(:,ja)
      end if
    end do
    if (ka > Faro_max_ex2a) stop 'exceeded max double excitations'
  end do

  call mma_deallocate(f)

end subroutine sigma2

subroutine sigma3(g,sgm,psi,ibsta,ibend)
! sigma3(ia,ib) = sum_ja,jb sum_tu,vx <jb|E_tu|ib> <ja|E_vx|ia> g_tuvx C(ja,jb)

  ! integrals
  real(kind=wp), intent(in) :: g(Faro_nOrb,Faro_nOrb,Faro_nOrb,Faro_nOrb)
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

  call mma_allocate(ja,Faro_max_LRs,label='ja')
  call mma_allocate(ia,Faro_max_LRs,label='ia')
  call mma_allocate(sgn_vx,Faro_max_LRs,label='sgn_vx')
  call mma_allocate(Ctmp,Faro_max_LRs,Faro_ndetb,label='Ctmp')
  call mma_allocate(Vtmp,Faro_max_LRs,label='Vtmp')
  call mma_allocate(f,Faro_ndetb,label='f')

  do v=1,Faro_nOrb
    do x=1,Faro_nOrb
      ! set up L(I), R(I), sgn(I) defined by L(ia) = E_tu R(ia)
      call LRs_init(v,x,Faro_nela,Faro_nOrb,ja,ia,sgn_vx,n_couples)
      do i=1,n_couples
        Ctmp(i,:) = psi(:,ja(i))*sgn_vx(i)
      end do
      do ib=ibsta,ibend
        f = Zero
        do tu=1,Faro_max_ex1b
          t = Faro_ex1_b(tu,ib)%p
          u = Faro_ex1_b(tu,ib)%q
          sgn_tu = Faro_ex1_b(tu,ib)%sgn
          jb = Faro_ex1_b(tu,ib)%rank
          ! for singlets, only tu >= vx are used
          if (Faro_mult == 1) then
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
        do tu=1,Faro_max_ex1b
          jb = Faro_ex1_b(tu,ib)%rank
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
        if (kb > Faro_max_ex1b) stop 'exceeded max single excitations'
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

subroutine Faro_ex1_init(k,n,ex1_table)

  use Faroald_SecQuant, only: Faro_ex1, Faro_fase, Faro_lex_init, Faro_lex_next, Faro_lexrank

  integer(kind=iwp), intent(in) :: k, n
  type(ex1_struct), intent(out) :: ex1_table(:,:)
  integer(kind=iwp) :: counter, det, idet, Faro_nDet, p, pq, q, tmp

  !write(u6,'(1x,a)') 'excitation table'
  !write(u6,'(1x,a)') 'p   q   I   J'
  Faro_nDet = Faro_Binom(k,n)
  det = Faro_lex_init(k,n)
  counter = 0
  do idet=1,Faro_nDet
    pq = 0
    do p=1,Faro_nOrb
      do q=1,Faro_nOrb
        tmp = Faro_ex1(p,q,det)
        if (tmp /= -1) then
          pq = pq+1
          ex1_table(pq,idet)%p = p
          ex1_table(pq,idet)%q = q
          ex1_table(pq,idet)%sgn = Faro_fase(tmp)
          ex1_table(pq,idet)%rank = Faro_lexrank(tmp)
          counter = counter+1
          !write(u6,'(1x,4i4)') p,q,idet,Faro_fase(tmp)*Faro_lexrank(tmp)
        end if
      end do
    end do
    det = Faro_lex_next(det)
  end do

end subroutine Faro_ex1_init

subroutine Faro_hDiag(h,g,diag,CSFDIA,nConf)

  use Constants, only: Zero

  integer(kind=iwp), intent(in):: nConf
  real(kind=wp) :: CSFDIA(nConf)
  real(kind=wp), intent(in) :: h(Faro_nOrb,Faro_nOrb), g(Faro_nOrb,Faro_nOrb,Faro_nOrb,Faro_nOrb)

  ! one element per spin combination (ib <= ia for singlets, see below)
  real(kind=wp), intent(out) :: diag(Faro_nComb_tot)

  real(kind=wp), allocatable :: Da(:), Db(:)

  integer(kind=iwp), allocatable :: occa(:,:), occb(:,:)

  integer(kind=iwp) :: ia, ib, i, j, k, iaLow
! ToDo cleanup this!
! integer(kind=iwp) :: MinOpen, MaxOpen, nOpen, nAEl, nBEl, iDet, jDet, iCSF, iCNF, jCNF, iCSFOff
  integer(kind=iwp) :: MinOpen, MaxOpen, nOpen, nAEl,       iDet, jDet, iCSF, iCNF, jCNF, iCSFOff
  integer(kind=iwp) :: SumDet, SumCSF, ipat, ipat_match
  real(kind=wp) :: EAver, EDet

  integer(kind=iwp) :: occ(Faro_nOrb)

  call mma_allocate(Da,Faro_ndeta,label='Da')
  call mma_allocate(Db,Faro_ndetb,label='Db')

  call mma_allocate(occa,Faro_nela,Faro_ndeta,label='OccA')
  call mma_allocate(occb,Faro_nelb,Faro_ndetb,label='OccB')

  call build_occ_alpha(occa)
  call build_occ_beta (occb)
  !
  ! Alpha-string energies
  !
  do ia=1,Faro_ndeta
     Da(ia) = Zero
     do i=1,Faro_nela
        Da(ia) = Da(ia) + h(occa(i,ia),occa(i,ia))
     end do
     do i=1,Faro_nela
        do j=i+1,Faro_nela
           Da(ia) = Da(ia)                                     &
              + g(occa(i,ia),occa(i,ia),occa(j,ia),occa(j,ia)) &
              - g(occa(i,ia),occa(j,ia),occa(j,ia),occa(i,ia))
        end do
     end do
  end do
  !
  ! Beta-string energies
  !
  do ib=1,Faro_ndetb
     Db(ib) = Zero
     do i=1,Faro_nelb
        Db(ib) = Db(ib) + h(occb(i,ib),occb(i,ib))
     end do
     do i=1,Faro_nelb
        do j=i+1,Faro_nelb
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
  do ib=1,Faro_ndetb
     iaLow = 1
     if (Faro_mult==1) iaLow=ib
     do ia=iaLow,Faro_ndeta
        k = k + 1
!       Write (6,*) 'k=',k,' -> ib, ia:',ib,ia
        diag(k) = Da(ia) + Db(ib)
        do i=1,Faro_nela
           do j=1,Faro_nelb
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

  MinOpen=Faro_mult-1
  MaxOpen=Min(Faro_nEl,Faro_nOrb)
  If (mod(MaxOpen-Faro_nEl,2) /= 0) MaxOpen=MaxOpen-1

! nOpen = Index of open-shell class. Equal the actual number of open shells.
  CSFDIA(:)=Zero
  iCSFOff=0
  sumDet = 0
  sumCSF = 0
  Do nOpen = MinOpen, MaxOpen, 2
     nAEl=(nOpen+(Faro_mult-1))/2
!    nBEl=(nOpen-(Faro_mult-1))/2
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
!    iDet = Faro_nComb(ipat_match)
!
! rather than the full determinant multiplicity
!
!    binom(nOpen,nAEl)
!
! when constructing the configuration average.
!
     iDet = Faro_Binom(nAEl,nOpen)
     iCSF=ncsf_open_shells(nOpen,Faro_mult-1)
     iCNF=nconf_open_shells(Faro_nOrb,Faro_nEl,nOpen)
!    Write(6,'(A,5I6)') 'nOpen, nAEl, nBEl, nDet, nCSF=', nOpen, nAEl, nBEl, iDet, iCSF
     sumDet = sumDet + iDet*iCNF
     sumCSF = sumCSF + iCSF*iCNF
     Do jCNF=1, iCNF

        call configuration_from_index(nOpen,jCNF,occ)

        ipat_match = 0

        do ipat=1,Faro_npat
           if (all(occ(:) == Faro_occpat(:,ipat))) then
              ipat_match = ipat
              exit
           end if
        end do
!       write(u6,'(A,3I6,A,20I2)') 'nOpen,jCNF,ipat=',nOpen,jCNF,ipat_match,' occ=',occ(:)

        EAver = Zero

        iDet = Faro_nComb(ipat_match)
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

   nOpen = Faro_nsing(ipat)

   if (nOpen == 0) then

      call Faro_combination_diagonal(ipat,h,g,e)

   else

      call mma_allocate(comb,nOpen,Faro_nComb(ipat),label='Comb')

      call Faro_pattern_combinations(ipat,nOpen,nCombPat,comb)

      call Faro_combination_diagonal(ipat,h,g,e,comb(:,jDet))

      call mma_deallocate(comb)

   end if

end subroutine


  subroutine build_occ_alpha(occ)

    integer(kind=iwp), intent(out) :: occ(Faro_nela,Faro_ndeta)

    integer(kind=iwp) :: ia, ipq, nocc

    do ia=1,Faro_ndeta

       nocc = 0

       do ipq=1,Faro_max_ex1a

          if (Faro_ex1_a(ipq,ia)%p == Faro_ex1_a(ipq,ia)%q .and. &
              Faro_ex1_a(ipq,ia)%rank == ia) then

             nocc = nocc + 1
             occ(nocc,ia) = Faro_ex1_a(ipq,ia)%p

          end if

       end do

       if (nocc /= Faro_nela) then
          write(u6,*) 'build_occ_alpha: wrong occupation count'
          call Abend()
       end if

    end do

  end subroutine build_occ_alpha

  subroutine build_occ_beta(occ)

    integer(kind=iwp), intent(out) :: occ(Faro_nelb,Faro_ndetb)

    integer(kind=iwp) :: ib, ipq, nocc

    do ib=1,Faro_ndetb

       nocc = 0

       do ipq=1,Faro_max_ex1b

          if (Faro_ex1_b(ipq,ib)%p == Faro_ex1_b(ipq,ib)%q .and. &
              Faro_ex1_b(ipq,ib)%rank == ib) then

             nocc = nocc + 1
             occ(nocc,ib) = Faro_ex1_b(ipq,ib)%p

          end if

       end do

       if (nocc /= Faro_nelb) then
          write(u6,*) 'build_occ_beta: wrong occupation count'
          call Abend()
       end if

    end do

  end subroutine build_occ_beta


end subroutine Faro_hDiag

subroutine Faro_build_pattern_diagonal_approx(h,g,diag)

! Occupation-pattern (Olsen-type) diagonal.
! Does NOT represent the exact CSF diagonal.

  use Constants, only: Zero

  real(kind=wp), intent(in) :: h(Faro_nOrb,Faro_nOrb), g(Faro_nOrb,Faro_nOrb,Faro_nOrb,Faro_nOrb)
  real(kind=wp), intent(out) :: diag(Faro_npat)

  integer(kind=iwp) :: ipat, p, q

  integer(kind=iwp) :: np, nq, exch

  do ipat=1,Faro_npat

     diag(ipat) = Zero

     !
     ! One-electron contribution
     !
     do p=1,Faro_nOrb

        np = Faro_occpat(p,ipat)

        diag(ipat) = diag(ipat) + real(np,wp)*h(p,p)

     end do

     !
     ! Two-electron contribution
     !
     do p=1,Faro_nOrb-1

        np = Faro_occpat(p,ipat)
        if (np == 0) cycle

        do q=p+1,Faro_nOrb

           nq = Faro_occpat(q,ipat)
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

end subroutine Faro_build_pattern_diagonal_approx

subroutine det_occ_pattern(ia,ib,occ)

  integer(kind=iwp), intent(in) :: ia, ib

  integer(kind=iwp), intent(out) :: occ(Faro_nOrb)

  integer(kind=iwp) :: p, pq

  occ(:) = 0

  !
  ! Alpha occupations
  !

  do pq=1,Faro_max_ex1a

     if (Faro_ex1_a(pq,ia)%p == Faro_ex1_a(pq,ia)%q .and. &
         Faro_ex1_a(pq,ia)%rank == ia) then

        p = Faro_ex1_a(pq,ia)%p
        occ(p) = occ(p) + 1

     end if

  end do

  !
  ! Beta occupations
  !

  do pq=1,Faro_max_ex1b

     if (Faro_ex1_b(pq,ib)%p == Faro_ex1_b(pq,ib)%q .and. &
         Faro_ex1_b(pq,ib)%rank == ib) then

        p = Faro_ex1_b(pq,ib)%p
        occ(p) = occ(p) + 1

     end if

  end do

end subroutine det_occ_pattern

subroutine Faro_verify_occ_patterns()

  integer(kind=iwp) :: ia, ib, pq, p

  integer(kind=iwp) :: occ(Faro_nOrb)

  integer(kind=iwp) :: na, nb

  do ib=1,Faro_ndetb

     do ia=1,Faro_ndeta

        occ(:) = 0

        !
        ! Recover alpha occupations
        !

        na = 0

        do pq=1,Faro_max_ex1a

           if (Faro_ex1_a(pq,ia)%p == Faro_ex1_a(pq,ia)%q .and. &
               Faro_ex1_a(pq,ia)%rank == ia) then

              p = Faro_ex1_a(pq,ia)%p

              occ(p) = occ(p) + 1
              na = na + 1

           end if

        end do

        !
        ! Recover beta occupations
        !

        nb = 0

        do pq=1,Faro_max_ex1b

           if (Faro_ex1_b(pq,ib)%p == Faro_ex1_b(pq,ib)%q .and. &
               Faro_ex1_b(pq,ib)%rank == ib) then

              p = Faro_ex1_b(pq,ib)%p

              occ(p) = occ(p) + 1
              nb = nb + 1

           end if

        end do

        !
        ! Consistency checks
        !

        if (na /= Faro_nela) then
           write(u6,*) 'verify_occ_patterns: wrong alpha count'
           write(u6,*) 'ia=',ia,' na=',na,' nela=',Faro_nela
           call Abend()
        end if

        if (nb /= Faro_nelb) then
           write(u6,*) 'verify_occ_patterns: wrong beta count'
           write(u6,*) 'ib=',ib,' nb=',nb,' nelb=',Faro_nelb
           call Abend()
        end if

        if (sum(occ) /= Faro_nEl) then
           write(u6,*) 'verify_occ_patterns: wrong electron count'
           write(u6,*) 'ia=',ia,' ib=',ib
           write(u6,*) 'sum(occ)=',sum(occ)
           write(u6,*) 'my_nel  =',Faro_nEl
           call Abend()
        end if

        if (maxval(occ) > 2) then
           write(u6,*) 'verify_occ_patterns: occupation > 2'
           write(u6,*) 'ia=',ia,' ib=',ib
           call Abend()
        end if

     end do

  end do

end subroutine Faro_verify_occ_patterns

subroutine Faro_build_patterns()

  use stdalloc, only: mma_allocate, mma_deallocate

  integer(kind=iwp) :: ia, ib, k, ipat, jpat

  integer(kind=iwp) :: occ(Faro_nOrb)

  integer(kind=iwp), allocatable :: tmp(:,:)

  logical :: found

  !
  ! Worst case: every determinant product generates
  ! its own occupation pattern.
  !

  call mma_allocate(Faro_occpat,Faro_nOrb,Faro_nDet,label='OccPat')
  call mma_allocate(Faro_ipat_of_det,Faro_nDet,label='IPatDet')

  Faro_npat = 0
  k    = 0

  do ib=1,Faro_ndetb

     do ia=1,Faro_ndeta

        k = k + 1

        call det_occ_pattern(ia,ib,occ)

        found = .false.

        do jpat=1,Faro_npat

           if (all(occ(:) == Faro_occpat(:,jpat))) then

              Faro_ipat_of_det(k) = jpat
              found = .true.

              exit

           end if

        end do

        if (.not. found) then

           Faro_npat = Faro_npat + 1

           Faro_occpat(:,Faro_npat) = occ(:)

           Faro_ipat_of_det(k) = Faro_npat

        end if

     end do

  end do

  !
  ! Compress Faro_occpat to its actual size.
  !

  call mma_allocate(tmp,Faro_nOrb,Faro_npat,label='TmpPat')

  tmp(:,:) = Faro_occpat(:,1:Faro_npat)

  call mma_deallocate(Faro_occpat)

  call mma_allocate(Faro_occpat,Faro_nOrb,Faro_npat,label='OccPat')

  Faro_occpat(:,:) = tmp(:,:)

  call mma_deallocate(tmp)

  !
  ! Consistency checks
  !

  do k=1,Faro_nDet
     ipat = Faro_ipat_of_det(k)

     if (ipat < 1 .or. ipat > Faro_npat) then
        write(u6,*) 'build_patterns: invalid pattern index'
        write(u6,*) 'k=',k,' ipat=',ipat,' npat=',Faro_npat
        call Abend()
     end if

  end do

#ifdef _DEBUGPRINT_
  do k = 1, Faro_npat
     write(u6,'(A,20I2)') 'occ=',Faro_occpat(:,k)
  end do
  write(u6,'(A,I10)') ' NDet = ', Faro_nDet
  write(u6,'(A,I10)') ' NPat = ', Faro_npat
  write(u6,'(A,F12.4)') ' Compression = ', real(Faro_nDet,wp)/real(Faro_npat,wp)
#endif

#ifdef _DEBUGPRINT_
do ipat=1,Faro_npat
   write(u6,'(A,I4,A,20I2)') &
        'ipat=',ipat, &
        ' occ=',Faro_occpat(:,ipat)
end do
#endif


end subroutine Faro_build_patterns

!#define _Print_Combinations_
subroutine Faro_analyse_patterns()

  integer(kind=iwp) :: ipat, p

#ifdef _Print_Combinations_
integer(kind=iwp) :: nopen, ncomb_pat, icomb
integer(kind=iwp), allocatable :: Comb(:,:)
integer(kind=iwp) :: occa_tmp(Faro_nela), occb_tmp(Faro_nelb)
#endif

  call mma_allocate(Faro_ndoub,Faro_npat,label='NDoub')
  call mma_allocate(Faro_nsing,Faro_npat,label='NSing')
  call mma_allocate(Faro_nSpin_Comb,Faro_npat,label='NSpinComb')
  call mma_allocate(Faro_nComb,Faro_npat,label='NComb')

  do ipat=1,Faro_npat

     Faro_ndoub(ipat) = 0
     Faro_nsing(ipat) = 0

     do p=1,Faro_nOrb

        select case (Faro_occpat(p,ipat))

        case (2)

           Faro_ndoub(ipat) = Faro_ndoub(ipat) + 1

        case (1)

           Faro_nsing(ipat) = Faro_nsing(ipat) + 1

        end select

     end do

     if (2*Faro_ndoub(ipat)+Faro_nsing(ipat) /= Faro_nEl) then

        write(u6,*) 'analyse_patterns: electron count error'
        write(u6,*) 'ipat=',ipat

        call Abend()

     end if
!
! Lucia uses spin combinations only for singlets (MS2 = 0).
! For all other spins the code falls back to the ordinary
! SD representation. Therefore Faro_nComb equals the number of
! alpha/beta assignments among the open shells.
!
     Faro_nSpin_Comb(ipat) = spin_comb_count(Faro_nsing(ipat),Faro_mult-1)
     if (Faro_mult /= 1) then
       !
       ! Lucia uses spin combinations only for singlets.
       ! For non-singlets the code falls back to the ordinary
       ! SD representation
       !
       Faro_nComb(ipat) = Faro_nSpin_Comb(ipat)
     else
       if (Faro_nsing(ipat) == 0) then
         Faro_nComb(ipat) = 1
       else
         Faro_nComb(ipat) = Faro_nSpin_Comb(ipat)/2
       end if
     end if

#ifdef _Print_Combinations_
   call mma_allocate(comb,Faro_nsing(ipat),Faro_nSpin_Comb(ipat),label='Comb')
   call Faro_pattern_combinations(ipat,nopen,ncomb_pat,comb)

   do icomb=1,ncomb_pat

      call Faro_combination_occupations(ipat,comb(:,icomb),occa_tmp,occb_tmp)

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

do ipat=1,Faro_npat

   write(u6,'(5I8)') &
        ipat, &
        Faro_nsing(ipat), &
        count(Faro_ipat_of_det == ipat), &
        Faro_nComb(ipat), &
        ncsf_open_shells(Faro_nsing(ipat),Faro_mult-1)

end do
#endif

call mma_allocate(Faro_ibComb,Faro_npat,label='IBComb')

Faro_ibComb(1) = 1

do ipat=2,Faro_npat
   Faro_ibComb(ipat) = Faro_ibComb(ipat-1) + Faro_nComb(ipat-1)
end do

if (Faro_mult == 1) then
   Faro_nComb_tot = Faro_ndeta*(Faro_ndeta+1)/2
else
   Faro_nComb_tot = Faro_ndeta*Faro_ndetb
end if

call Faro_build_ictsdt()


end subroutine Faro_analyse_patterns

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
     spin_comb_count = Faro_Binom(iael,iopen)
  end if

end function spin_comb_count

subroutine spncom_faroald(nopen,ms2,Faro_nComb,comb)

  integer(kind=iwp), intent(in) :: nopen, ms2
  integer(kind=iwp), intent(in) :: Faro_nComb
  integer(kind=iwp), intent(out) ::  comb(nopen,Faro_nComb)
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
       if (Faro_mult /= 1 .or. work(1) == 1) then
         icomb = icomb + 1
         comb(:,icomb) = work(:)
       end if
     end if

  end do

end subroutine spncom_faroald


subroutine Faro_pattern_combinations(ipat,nopen,ncomb_pat,comb)

  integer(kind=iwp), intent(in) :: ipat
  integer(kind=iwp), intent(out) :: nopen, ncomb_pat
  integer(kind=iwp), intent(out) :: comb(:,:)

  integer(kind=iwp) :: p

  nopen = 0

  do p=1,Faro_nOrb
     if (Faro_occpat(p,ipat) == 1) nopen = nopen + 1
  end do

  ncomb_pat = Faro_nComb(ipat)

  call spncom_faroald(nopen,Faro_mult-1,ncomb_pat,comb)

end subroutine Faro_pattern_combinations

subroutine Faro_combination_occupations(ipat,comb,occa,occb)

  integer(kind=iwp), intent(in) :: ipat
  integer(kind=iwp), intent(in), optional :: comb(:)

  integer(kind=iwp), intent(out) :: occa(Faro_nela)
  integer(kind=iwp), intent(out) :: occb(Faro_nelb)

  integer(kind=iwp) :: ia
  integer(kind=iwp) :: ib
  integer(kind=iwp) :: ic
  integer(kind=iwp) :: p

  ia = 0
  ib = 0
  ic = 0

  do p=1,Faro_nOrb

     select case (Faro_occpat(p,ipat))

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

end subroutine Faro_combination_occupations

subroutine Faro_combination_diagonal(ipat,h,g,e,comb)

  real(kind=wp), intent(in) :: h(Faro_nOrb,Faro_nOrb)
  real(kind=wp), intent(in) :: g(Faro_nOrb,Faro_nOrb,Faro_nOrb,Faro_nOrb)

  integer(kind=iwp), intent(in) :: ipat
  integer(kind=iwp), intent(in), optional :: comb(:)

  real(kind=wp), intent(out) :: e

  integer(kind=iwp) :: occa(Faro_nela)
  integer(kind=iwp) :: occb(Faro_nelb)

  integer(kind=iwp) :: i,j

  call Faro_combination_occupations(ipat,comb,occa,occb)

  e = Zero

  do i=1,Faro_nela
     e = e + h(occa(i),occa(i))
  end do

  do i=1,Faro_nelb
     e = e + h(occb(i),occb(i))
  end do

  do i=1,Faro_nela
     do j=i+1,Faro_nela
        e = e + g(occa(i),occa(i),occa(j),occa(j))
        e = e - g(occa(i),occa(j),occa(j),occa(i))
     end do
  end do

  do i=1,Faro_nelb
     do j=i+1,Faro_nelb
        e = e + g(occb(i),occb(i),occb(j),occb(j))
        e = e - g(occb(i),occb(j),occb(j),occb(i))
     end do
  end do

  do i=1,Faro_nela
     do j=1,Faro_nelb
        e = e + g(occa(i),occa(i),occb(j),occb(j))
     end do
  end do

end subroutine Faro_combination_diagonal

subroutine Faro_build_pattern_diagonal(h,g,diag)

  use Constants, only: Zero

  real(kind=wp), intent(in) :: h(Faro_nOrb,Faro_nOrb), g(Faro_nOrb,Faro_nOrb,Faro_nOrb,Faro_nOrb)
  real(kind=wp), intent(out) :: diag(Faro_npat)

  integer(kind=iwp) :: ipat, p,q, np,nq, exch


  do ipat=1,Faro_npat
     diag(ipat) = Zero

     !
     ! One-electron contribution
     !

     do p=1,Faro_nOrb
        np = Faro_occpat(p,ipat)
        diag(ipat) = diag(ipat) + real(np,wp)*h(p,p)
     end do

     !
     ! Two-electron contribution
     !

     do p=1,Faro_nOrb-1
        np = Faro_occpat(p,ipat)
        if (np == 0) cycle
        do q=p+1,Faro_nOrb
           nq = Faro_occpat(q,ipat)
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

     do p=1,Faro_nOrb
        if (Faro_occpat(p,ipat) == 2) then
           diag(ipat) = diag(ipat) + g(p,p,p,p)
        end if
     end do
  end do

end subroutine Faro_build_pattern_diagonal

subroutine Faro_build_ictsdt()

  integer(kind=iwp) :: ia
  integer(kind=iwp) :: ib
  integer(kind=iwp) :: iaLow
  integer(kind=iwp) :: k

  call mma_allocate(Faro_ictsdt,Faro_nComb_tot,label='ICTSDT')

  k = 0

  do ib=1,Faro_ndetb

     iaLow = 1
     if (Faro_mult == 1) iaLow = ib

     do ia=iaLow,Faro_ndeta

        k = k + 1

        Faro_ictsdt(k) = k

     end do

  end do

end subroutine Faro_build_ictsdt


integer(kind=iwp) function nconf_open_shells(norb,nel,nopen)

  integer(kind=iwp), intent(in) :: norb
  integer(kind=iwp), intent(in) :: nel
  integer(kind=iwp), intent(in) :: nopen

  integer(kind=iwp) :: Faro_ndoub

  if (mod(nel-nopen,2) /= 0) then
     nconf_open_shells = 0
     return
  end if

  Faro_ndoub = (nel-nopen)/2

  nconf_open_shells = &
       Faro_Binom(nopen,norb) * &
       Faro_Binom(Faro_ndoub,norb-nopen)

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

  ncsf_open_shells = Faro_Binom(k,nopen)

  if (k > 0) then
     ncsf_open_shells = ncsf_open_shells - &
                        Faro_Binom(k-1,nopen)
  end if

end function ncsf_open_shells


subroutine configuration_from_index(nOpen,jCNF,occ)

  integer(kind=iwp), intent(in)  :: nOpen
  integer(kind=iwp), intent(in)  :: jCNF
  integer(kind=iwp), intent(out) :: occ(Faro_nOrb)

  integer(kind=iwp) :: iconf_exp(Faro_nEl)
  integer(kind=iwp) :: occ_local(Faro_nOrb)

  integer(kind=iwp) :: i
  integer(kind=iwp) :: nopen_local
  integer(kind=iwp) :: count_conf
  integer(kind=iwp) :: nonew

  occ(:) = 0

  !
  ! Initial configuration
  !
  call nxt_conf_faroald(iconf_exp,Faro_nEl,Faro_nOrb,1,nonew)

  count_conf = 0

  do while (nonew == 0)

     !
     ! Expanded occupation list -> occupation pattern
     !
     occ_local(:) = 0

     do i=1,Faro_nEl
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
     call nxt_conf_faroald(iconf_exp,Faro_nEl,Faro_nOrb,0,nonew)

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

subroutine LRs_init(p,q,Faro_nEl,Faro_nOrb,L,R,sgn,counter)
! for a pair of orbitals p and q, and determinants
! generated by Faro_nEl electrons in Faro_nOrb spin orbitals,
! re-enumerate all non-vanishing couples connected
! through: jdet = E_pq idet. The number of couples
! is given by n_det, L(i) = jdet and R(i) = idet,
! where i is a counter from 1 to n_det. The sign of
! jdet is stored in the sgn array.

  use Faroald_SecQuant, only: Faro_ex1, Faro_fase, Faro_lex_init, Faro_lex_next, Faro_lexrank

  integer(kind=iwp), intent(in) :: p, q, Faro_nEl, Faro_nOrb
  integer(kind=iwp), intent(out) :: L(:), R(:), sgn(:), counter
  integer(kind=iwp) :: det, idet, Faro_nDet, tmp

  Faro_nDet = Faro_Binom(Faro_nEl,Faro_nOrb)
  det = Faro_lex_init(Faro_nEl,Faro_nOrb)
  counter = 0
  do idet=1,Faro_nDet
    tmp = Faro_ex1(p,q,det)
    if (tmp /= -1) then
      counter = counter+1
      R(counter) = idet
      sgn(counter) = Faro_fase(tmp)
      L(counter) = Faro_lexrank(tmp)
    end if
    det = Faro_lex_next(det)
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

end module Faroald_Core