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
! Copyright (C) 2022, Arta Safari                                      *
!***********************************************************************

module spin_correlation

use CI_solver_util, only: rdm_from_runfile
use rasscf_global, only: iAdr15, lRoots, nacpar, nacpr2
use index_symmetry, only: one_el_idx_flatten, two_el_idx_flatten
use general_data, only: JobIPH
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero
use Definitions, only: wp, iwp, u6

implicit none
private

integer(kind=iwp) :: same_orbs
logical(kind=iwp) :: tRootGrad = .false.
integer(kind=iwp), allocatable :: orb_range_p(:), orb_range_q(:)

public :: orb_range_p, orb_range_q, same_orbs, spin_correlation_driver, tRootGrad

contains

subroutine spin_correlation_driver(orb_range_p,orb_range_q,iroot)
!! spin-spin-correlation function using orbital-resolved 2RDMs.
!! For details see Dobrautz et al. 2021, 10.1021/acs.jctc.1c00589.

  integer(kind=iwp), intent(in) :: orb_range_p(:), orb_range_q(:), iroot(:)
  integer(kind=iwp) :: i, j, jDisk
  logical(kind=iwp) :: disk_pointer_moved
  real(kind=wp), allocatable :: dmat(:), dspn(:), pamat(:), psmat(:), spin_correlations(:)

  jDisk = iAdr15(3)
  call mma_allocate(dmat,nacpar,Label='dmat')
  call mma_allocate(dspn,nacpar,Label='dspn')
  call mma_allocate(pamat,nacpr2,Label='pamat')
  call mma_allocate(psmat,nacpr2,Label='psmat')
  call mma_allocate(spin_correlations,size(iroot))
  spin_correlations(:) = Zero

  ! the following disk pointer logic is necessary, because a user may
  ! optimise less roots than included in the Davidson matrix, e.g.
  ! roots (1, 2, 5, 8, 10). In &RASSCF RDMs for all Lroots
  ! will be computed, hence the memory reference here needs to be
  ! incremented accordingly to compute the SSCR for the correct roots.

  write(u6,'(a)') new_line('a')
  do i=1,LRoots
    disk_pointer_moved = .false.
    do j=1,size(iroot)
      if (iroot(j) == i) then
        call rdm_from_runfile(dmat,dspn,psmat,pamat,jDisk)
        disk_pointer_moved = .true.
        spin_correlations(j) = correlation_func(orb_range_p,orb_range_q,dmat,psmat,pamat)
        write(u6,'(a,i2,a,f12.8)') '::    RASSCF root number ',iroot(j),' Spin Correlation:  ',spin_correlations(j)
      end if
    end do
    if (.not. disk_pointer_moved) then
      ! dummy write to move jDisk pointer
      call ddafile(JOBIPH,0,dmat,NACPAR,jDisk)
      call ddafile(JOBIPH,0,dspn,NACPAR,jDisk)
      call ddafile(JOBIPH,0,psmat,NACPR2,jDisk)
      call ddafile(JOBIPH,0,pamat,NACPR2,jDisk)
    end if
  end do

  ! for testing purposes
  call Add_Info('spin correlation',spin_correlations(1),1,8)

  call mma_deallocate(dmat)
  call mma_deallocate(dspn)
  call mma_deallocate(pamat)
  call mma_deallocate(psmat)
  call mma_deallocate(spin_correlations)

end subroutine spin_correlation_driver

function correlation_func(orb_range_p,orb_range_q,dmat,psmat,pamat) result(corr)
!! extract spin-spin-correlation function from orbital resolved RDMs.

  use Constants, only: Half

  real(kind=wp) corr
  integer(kind=iwp), intent(in) :: orb_range_p(:), orb_range_q(:)
  real(kind=wp), intent(in) :: dmat(nacpar), psmat(nacpr2), pamat(nacpr2)
  integer(kind=iwp) :: p, pp, pppp, ppqq, pqqp, q, rp, rq
  real(kind=wp) :: onerdm_pp, twordm_pppp, twordm_ppqq, twordm_pqqp

  corr = Zero

  do p=1,size(orb_range_p)
    do q=1,size(orb_range_q)
      ! dummy variables to save space
      rp = orb_range_p(p)
      rq = orb_range_q(q)
      if (rp /= rq) then
        pqqp = two_el_idx_flatten(rp,rq,rq,rp)
        ppqq = two_el_idx_flatten(rp,rp,rq,rq)
        twordm_pqqp = psmat(pqqp)-pamat(pqqp)
        twordm_ppqq = 2*(psmat(ppqq)+pamat(ppqq))

        corr = corr-Half*(twordm_pqqp+Half*twordm_ppqq)
      else
        pppp = two_el_idx_flatten(rp,rp,rp,rp)
        pp = one_el_idx_flatten(rp,rp)

        twordm_pppp = 2*(psmat(pppp)+pamat(pppp))
        onerdm_pp = dmat(pp)
        corr = corr+0.75_wp*(onerdm_pp-twordm_pppp)
      end if
    end do
  end do

end function correlation_func

end module spin_correlation
