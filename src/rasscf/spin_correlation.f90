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
  use definitions, only: wp, u6
  use stdalloc, only: mma_allocate, mma_deallocate
  use CI_solver_util, only: rdm_from_runfile
  use rasscf_data, only : NRoots, iAdr15, nacpar, nacpr2
  use index_symmetry, only : one_el_idx_flatten, two_el_idx_flatten

  implicit none
  private
  public :: spin_correlation_driver
  integer, allocatable, public, save :: orb_range_p(:), orb_range_q(:)
  integer, public :: same_orbs


contains


  subroutine spin_correlation_driver(orb_range_p, orb_range_q, iroot)
    !! spin-spin-correlation function using orbital-resolved 2RDMs.
    !! For details see Dobrautz et al. 2021, 10.1021/acs.jctc.1c00589.
    integer, intent(in) :: orb_range_p(:), orb_range_q(:), iroot(:)
    real(wp), allocatable :: spin_correlations(:)
    real(wp) :: dmat(nacpar), dspn(nacpar), pamat(nacpr2), psmat(nacpr2)
    integer :: jDisk, i, j

    jDisk = iAdr15(3)
    call mma_allocate(spin_correlations, size(iroot))
    spin_correlations(:) = 0.0_wp

    write(u6,'(a)') new_line('a')
    do i = 1, NRoots
      do j = 1, size(iroot)
        if (iroot(j) == i) then
          call rdm_from_runfile(dmat, dspn, psmat, pamat, jDisk)
          spin_correlations(j) = correlation_func(orb_range_p, orb_range_q, &
                                                  dmat, psmat, pamat)
          write(u6,'(a,i2,a,f12.8)') '::    RASSCF root number ', iroot(j), &
            ' Spin Correlation:  ', spin_correlations(j)
        end if
      end do
    end do

    ! for testing purposes
    call Add_Info('spin correlation', spin_correlations(1), 1, 8)

    call mma_deallocate(spin_correlations)
  end subroutine spin_correlation_driver


  real(wp) function correlation_func(orb_range_p, orb_range_q, &
                                     dmat, psmat, pamat) result(corr)
    !! extract spin-spin-correlation function from orbital resolved RDMs.
    integer, intent(in) :: orb_range_p(:), orb_range_q(:)
    real(wp), intent(in) :: dmat(nacpar), psmat(nacpr2), pamat(nacpr2)
    integer :: rp, rq, p, q, pp, pppp, pqqp, ppqq
    real(wp) :: twordm_pqqp, twordm_ppqq, twordm_pppp, onerdm_pp

    corr = 0.0_wp

    do p = 1, size(orb_range_p)
      do q = 1, size(orb_range_q)
        ! dummy variables to save space
        rp = orb_range_p(p); rq = orb_range_q(q)
        if (rp /= rq) then
          pqqp = two_el_idx_flatten(rp, rq, rq, rp)
          ppqq = two_el_idx_flatten(rp, rp, rq, rq)
          twordm_pqqp = psmat(pqqp) - pamat(pqqp)
          twordm_ppqq = 2 * (psmat(ppqq) + pamat(ppqq))

          corr = corr - 0.5_wp * (twordm_pqqp + 0.5_wp * twordm_ppqq)
        else
          pppp = two_el_idx_flatten(rp, rp, rp, rp)
          pp = one_el_idx_flatten(rp, rp)

          twordm_pppp = 2 * (psmat(pppp) + pamat(pppp))
          onerdm_pp = dmat(pp)
          corr = corr + 0.75_wp * (onerdm_pp - twordm_pppp)
        end if
      end do
    end do
  end function correlation_func

end module spin_correlation
