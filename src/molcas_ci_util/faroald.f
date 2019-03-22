************************************************************************
* This file is part of OpenMolcas.                                     *
*                                                                      *
* OpenMolcas is free software; you can redistribute it and/or modify   *
* it under the terms of the GNU Lesser General Public License, v. 2.1. *
* OpenMolcas is distributed in the hope that it will be useful, but it *
* is provided "as is" and without any express or implied warranties.   *
* For more details see the full text of the license in the file        *
* LICENSE or in <http://www.gnu.org/licenses/>.                        *
*                                                                      *
* Copyright (C) 2014, Steven Vancoillie                                *
************************************************************************
      module faroald
!     written by Steven Vancoillie, summer 2014
!
!     The faroald module handles sigma updates as the result of acting
!     with the hamiltonian operator on a CI vector: s = H c, where s and c
!     are CI expansions in determinant basis.
!
!     The implementation follows the minimum operation count algorithm
!     published by Olsen & Co in J. Chem. Phys. 89, 2185 (1988).

      implicit none
      save

      ! wavefunction info
      integer :: my_nel, my_norb, mult
      integer :: nela, nelb, nhoa, nhob
      integer :: ndeta, ndetb, my_ndet

      ! excitation tables
      type ex1_struct
        integer :: p, q
        integer :: sgn
        integer :: rank
      end type

      type(ex1_struct), allocatable :: ex1_a(:,:), ex1_b(:,:)
      integer :: max_ex1a, max_ex1b
      integer :: max_ex2a, max_ex2b
      integer :: max_LRs

#ifdef _PROF_
      integer*8 :: nflop
#endif

      contains

      subroutine sigma_update(h,g,sgm,psi)
      ! The sigma update routine performs the following operation:
      ! |sgm> = H |psi>, with the hamiltonian defined by
      ! H = sum_tu h(t,u) E_tu + sum_tuvx g(t,u,v,x) E_tuvx.

      use second_quantization
      implicit none

      ! integrals
      real*8, intent(in)  :: h(my_norb,my_norb)
      real*8, intent(in)  :: g(my_norb,my_norb,my_norb,my_norb)
      real*8, allocatable :: k(:,:)

      ! wavefunctions
      real*8, intent(in)    :: psi(:,:)
      real*8, intent(inout) :: sgm(:,:)

      real*8, allocatable :: psiT(:,:)
      real*8, allocatable :: sgmT(:,:)

      ! orbital indices
      integer :: t, u, v

      ! determinant index ranges
      integer :: iasta, iaend, ibsta, ibend

      integer :: ierr

      ! profiling
#ifdef _PROF_
      real*8 :: t1_cpu, t2_cpu, tot_cpu
      real*8 :: t1_wall, t2_wall, tot_wall
      real*8 :: walltime, flops
#endif

      ! Distributes a dimension over processes
      ! and returns the range of this process.
      call par_range(ndeta,iasta,iaend)
      call par_range(ndetb,ibsta,ibend)

#ifdef _PROF_
      ! initialize flop count
      nflop = 0

      call timing (t1_cpu,tot_cpu,t1_wall,tot_wall)
#endif

      ! Set the sigma vector to 0
      sgm = 0.0d0

      ! First, construct a new effective one-electron integral matrix:
      ! k_tu = h_tu - 1/2 sum_v g_tvvu, to be used with sigma1/sigma2.
      allocate (k(my_norb, my_norb))
      do u=1,my_norb
        do t=1,my_norb
          ! g_tvvu = g_vtvu
          k(t,u) = 0.0d0
          do v=1,my_norb
            k(t,u) = k(t,u) + g(v,t,v,u)
          end do
          k(t,u) = h(t,u) - 0.5d0 * k(t,u)
        end do
      end do

!     Second, for sigma2 and sigma3, we are better off with the transpose
!     of sgm and/or psi, so allocate and assign them here.
      allocate(psiT(ndetb,ndeta), stat=ierr)
      if (ierr /= 0) stop 'sigma_update: could not allocate psiT'
      call dtrans(ndeta,ndetb,psi,ndeta,psiT,ndetb)

!     Now the actual contributions to sigma are computed. For a singlet
!     (mult = 1), sigma2 is not computed and sigma3 will only do half the
!     work.  But at the end, the transpose of sigma has to be added to
!     sigma. This should be more efficient than computing everything, but
!     note that the transpose operation is also included in the timings
!     used to compute the flop efficiency.

      call sigma1(k,g,sgm,psi,ibsta,ibend)

      if (mult.ne.1) then
        ! we need efficient access to sgm by using the transpose
        allocate(sgmT(ndetb,ndeta), stat=ierr)
        if (ierr /= 0) stop 'sigma_update: could not allocate sgmT'
        call dtrans(ndeta,ndetb,sgm,ndeta,sgmT,ndetb)

        call sigma2(k,g,sgmT,psiT,iasta,iaend)

        call dtrans(ndetb,ndeta,sgmT,ndetb,sgm,ndeta)
        deallocate(sgmT)
      end if

      call sigma3(g,sgm,psiT,ibsta,ibend)

      ! sum over all processes
      call gadsum(sgm,ndeta*ndetb)

      if (mult.eq.1) then
        ! for Ms = 0 (only used for singlet), sgm := sgm + sgm^T
        call transadd(ndeta,sgm,ndeta)
      end if

      deallocate(psiT)
      deallocate(k)

#ifdef _PROF_
      call timing (t2_cpu,tot_cpu,t2_wall,tot_wall)

      walltime = t2_wall - t1_wall

      if (walltime /= 0.0d0) then
        flops = nflop / walltime
        write(6,'(1x,a,2(f10.3,a))') &
          & 'sigma update: ', &
          & walltime , ' s, ', &
          & flops * 1.0d-9 , ' Gflops.'
      end if
#endif

      end subroutine

      subroutine sigma1(k,g,sgm,psi,ibsta,ibend)
!     sigma1 = sum_jb sum_tu <jb|E_tu|ib> (h_kl - 1/2 sum_v <tv|vx>) C(ia,jb)
!        + 1/2 sum_jb sum_tuvx <jb|E_tu E_vx|ib> g_tuvx C(ia,jb)

      implicit none

      ! integrals
      real*8, intent(in)  :: k(my_norb,my_norb)
      real*8, intent(in)  :: g(my_norb,my_norb,my_norb,my_norb)

      ! wavefunctions
      real*8, intent(in)    :: psi(:,:)
      real*8, intent(inout) :: sgm(:,:)

      ! local variables
      integer :: ib, jb, kb
      integer :: ibsta, ibend
      integer :: t, u, v, x, tu, vx
      integer :: sgn_tu, sgn_vx

      integer :: ierr

      real*8, allocatable :: f(:)

      allocate(f(ndetb), stat=ierr)
      if (ierr /= 0) stop 'could not allocate f'

      do ib=ibsta,ibend
        ! f array construction
        f=0.0d0
        do tu=1,max_ex1b
          t      = ex1_b(tu,ib) % p
          u      = ex1_b(tu,ib) % q
          sgn_tu = ex1_b(tu,ib) % sgn
          kb     = ex1_b(tu,ib) % rank
          f(kb) = f(kb) + sgn_tu * k(t,u)
          do vx=1,max_ex1b
            v      = ex1_b(vx,kb) % p
            x      = ex1_b(vx,kb) % q
            sgn_vx = ex1_b(vx,kb) % sgn
            jb     = ex1_b(vx,kb) % rank
            f(jb) = f(jb) + 0.5d0 * sgn_tu * sgn_vx * g(v,x,t,u)
          end do
        end do
        ! sigma addition
        kb = 0
        do jb=1,ndetb
          if (f(jb) /= 0.0d0) then
            kb = kb + 1
#ifdef _PROF_
            nflop = nflop + 2 * ndeta
#endif
            call daxpy_(ndeta,f(jb),psi(1,jb),1,sgm(1,ib),1)
          end if
        end do
        if (kb > max_ex2b) stop 'exceeded max double excitations'
      end do

      deallocate(f)

      end subroutine

      subroutine sigma2(k,g,sgm,psi,iasta,iaend)
!     sigma2 = sum_ja sum_tu <ja|E_tu|ia> (h_kl - 1/2 sum_v <tv|vx>) C(ja,ib)
!        + 1/2 sum_ja sum_tuvx <ja|E_tu E_vx|ia> g_tuvx C(ja,ib)

      implicit none

      ! integrals
      real*8, intent(in)  :: k(my_norb,my_norb)
      real*8, intent(in)  :: g(my_norb,my_norb,my_norb,my_norb)

      ! wavefunctions
      real*8, intent(in)    :: psi(:,:)
      real*8, intent(inout) :: sgm(:,:)

      ! local variables
      integer :: ia, ja, ka
      integer :: iasta, iaend
      integer :: t, u, v, x, tu, vx
      integer :: sgn_tu, sgn_vx

      integer :: ierr

      real*8, allocatable :: f(:)

      allocate(f(ndeta), stat=ierr)
      if (ierr /= 0) stop 'could not allocate f'

      do ia=iasta,iaend
        ! f array construction
        f=0.0d0
        do tu=1,max_ex1a
          t      = ex1_a(tu,ia) % p
          u      = ex1_a(tu,ia) % q
          sgn_tu = ex1_a(tu,ia) % sgn
          ka     = ex1_a(tu,ia) % rank
          f(ka) = f(ka) + sgn_tu * k(t,u)
          do vx=1,max_ex1a
            v      = ex1_a(vx,ka) % p
            x      = ex1_a(vx,ka) % q
            sgn_vx = ex1_a(vx,ka) % sgn
            ja     = ex1_a(vx,ka) % rank
            f(ja) = f(ja) + 0.5d0 * sgn_tu * sgn_vx * g(v,x,t,u)
          end do
        end do
        ! sigma addition
        ka = 0
        do ja=1,ndeta
          if (f(ja) /= 0.0d0) then
            ka = ka + 1
#ifdef _PROF_
            nflop = nflop + 2 * ndeta
#endif
            call daxpy_(ndetb,f(ja),psi(1,ja),1,sgm(1,ia),1)
          end if
        end do
        if (ka > max_ex2a) stop 'exceeded max double excitations'
      end do

      deallocate(f)

      end subroutine

      subroutine sigma3(g,sgm,psi,ibsta,ibend)
!     sigma3(ia,ib) = sum_ja,jb sum_tu,vx <jb|E_tu|ib> <ja|E_vx|ia> g_tuvx C(ja,jb)

      implicit none

      ! integrals
      real*8, intent(in)  :: g(my_norb,my_norb,my_norb,my_norb)

      ! wavefunctions
      real*8, intent(in)    :: psi(:,:)
      real*8, intent(inout) :: sgm(:,:)

      ! determinant indices
      integer, allocatable :: ia(:), ja(:)
      integer :: i, n_couples
      integer :: ib, jb, kb
      integer :: ibsta, ibend

      ! orbital indices
      integer :: t, u, v, x, tu
      integer, allocatable :: sgn_vx(:)
      integer :: sgn_tu

      real*8, allocatable :: f(:)
      real*8, allocatable :: Ctmp(:,:), Vtmp(:)

      integer :: ierr

      allocate(ja(max_LRs), ia(max_LRs), sgn_vx(max_LRs), stat=ierr)
      if (ierr /= 0) stop 'could not allocate L/R/sgn'

      allocate(Ctmp(max_LRs,ndetb), Vtmp(max_LRs), stat=ierr)
      if (ierr /= 0) stop 'could not allocate Ctmp/Vtmp'

      allocate(f(ndetb), stat=ierr)
      if (ierr /= 0) stop 'could not allocate f'

      do v=1,my_norb
        do x=1,my_norb
          ! set up L(I), R(I), sgn(I) defined by L(ia) = E_tu R(ia)
          call LRs_init(v, x, nela, my_norb, ja, ia, sgn_vx, n_couples)
          do i=1,n_couples
            Ctmp(i,:) = psi(:,ja(i)) * sgn_vx(i)
          end do
          do ib=ibsta,ibend
            f = 0.0d0
            do tu=1,max_ex1b
              t      = ex1_b(tu,ib) % p
              u      = ex1_b(tu,ib) % q
              sgn_tu = ex1_b(tu,ib) % sgn
              jb     = ex1_b(tu,ib) % rank
              ! for singlets, only tu >= vx are used
              if (mult.eq.1) then
                if (t < v) cycle
                if (t == v .and. u < x) cycle
                if (t == v .and. u == x) then
                  f(jb) = f(jb) + sgn_tu * 0.5d0 * g(t,u,v,x)
                  cycle
                end if
              end if
              f(jb) = f(jb) + sgn_tu * g(t,u,v,x)
            end do
            ! V(ia) = sum_jb f(jb) C'(ia,jb) for all ia
            Vtmp=0.0d0
            kb = 0
            ! loop over non-identical excitations
            do tu=1,max_ex1b
              jb = ex1_b(tu,ib) % rank
              if (jb /= ib .and. f(jb) /= 0.0d0) then
                kb = kb + 1
#ifdef _PROF_
                nflop = nflop + 2 * n_couples
#endif
                call daxpy_(n_couples,f(jb),Ctmp(1,jb),1,Vtmp,1)
              end if
            end do
            ! contribution from the identical excitations
            if (f(ib) /= 0.0d0) then
#ifdef _PROF_
              nflop = nflop + 2 * n_couples
#endif
              call daxpy_(n_couples,f(ib),Ctmp(1,ib),1,Vtmp,1)
            end if
            if (kb > max_ex1b) stop 'exceeded max single excitations'
            ! s3(R_ia,ib) = s3(R_ia,ib) + V(ia)
            do i=1,n_couples
              sgm(ia(i),ib) = sgm(ia(i),ib) + Vtmp(i)
            end do
          end do
        end do
      end do

      end subroutine

      subroutine ex1_init(k, n, ex1_table)
      use second_quantization
      implicit none
      integer, intent(in) :: k, n
      type(ex1_struct), intent(out) :: ex1_table(:,:)

      integer :: my_ndet, idet
      integer :: det, tmp
      integer :: p, q, pq

      integer :: counter

      !write(6,'(1x,a)') 'excitation table'
      !write(6,'(1x,a)') 'p   q   I   J'
      my_ndet = binom_coef(k,n)
      det = lex_init(k,n)
      counter = 0
      do idet=1,my_ndet
        pq = 0
        do p=1,my_norb
          do q=1,my_norb
            tmp = ex1(p,q,det)
            if (tmp.ne.-1) then
              pq = pq + 1
              ex1_table(pq,idet) % p    = p
              ex1_table(pq,idet) % q    = q
              ex1_table(pq,idet) % sgn  = fase(tmp)
              ex1_table(pq,idet) % rank = lexrank(tmp)
              counter = counter + 1
              !write(6,'(1x,4i4)') p, q, idet, fase(tmp)*lexrank(tmp)
            end if
          end do
        end do
        det=lex_next(det)
      end do
      end subroutine

      subroutine LRs_init(p, q, my_nel, my_norb, L, R, sgn, counter)
      ! for a pair of orbitals p and q, and determinants
      ! generated by my_nel electrons in my_norb spin orbitals,
      ! re-enumerate all non-vanishing couples connected
      ! through: jdet = E_pq idet. The number of couples
      ! is given by n_det, L(i) = jdet and R(i) = idet,
      ! where i is a counter from 1 to n_det. The sign of
      ! jdet is stored in the sgn array.
      use second_quantization
      implicit none
      integer, intent(in)  :: p, q, my_nel, my_norb
      integer, intent(out) :: L(:), R(:), sgn(:), counter

      integer :: my_ndet, idet
      integer :: det, tmp

      my_ndet = binom_coef(my_nel,my_norb)
      det = lex_init(my_nel,my_norb)
      counter = 0
      do idet=1,my_ndet
        tmp = ex1(p,q,det)
        if (tmp.ne.-1) then
          counter = counter + 1
          R(counter)   = idet
          sgn(counter) = fase(tmp)
          L(counter)   = lexrank(tmp)
        end if
        det=lex_next(det)
      end do
      end subroutine

      end module
