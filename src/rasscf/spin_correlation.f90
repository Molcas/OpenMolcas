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
* Copyright (C) 2022, Arta Safari                                      *
************************************************************************

module spin_correlation

  use definitions, only: wp
  use stdalloc, only: mma_allocate, mma_deallocate
  use general_data, only : nActEl
  use index_symmetry, only : two_el_idx

  implicit none

  private
  public :: spin_spin_correlation

contains

  pure function spin_spin_correlation(spinfree_2rdm, spinfree_1rdm, &
        orb_range_i, orb_range_j) result(spin_correlation)

    !! spin-spin-correlation function using orbital-resolved 2RDMs.
    !! For details see Dobrautz et al. 2021, 10.1021/acs.jctc.1c00589.

    real(dp), intent(in) :: spinfree_2rdm(:,:,:,:), spinfree_1rdm(:,:)
    integer, intent(in) :: orb_range_i(:), orb_range_j(:)

    real(dp) :: spin_correlation
    integer :: i, j

    spin_correlation = 0.0_dp

    ! spatial orbital labels p,q,r,s,...
    associate(ps => orb_range_i, qs => orb_range_j)
    do j = 1, size(ps)
      do i = 1, size(qs)
        if (ps(i) /= qs(j)) then
          spin_correlation = spin_correlation &
            - 0.5_dp * (spinfree_2rdm(ps(i),qs(j),qs(j),ps(i)) &
            + 0.5_dp * spinfree_2rdm(ps(i),ps(i),qs(j),qs(j)))
        else
          spin_correlation = spin_correlation &
            + 0.75_dp * (spinfree_1rdm(ps(i),ps(i)) &
            - spinfree_2rdm(ps(i),ps(i),ps(i),ps(i))) &
            - 0.5_dp * (spinfree_2rdm(ps(i),ps(j),ps(j),ps(i)) &
            + 0.5_dp * spinfree_2rdm(ps(i),ps(i),ps(j),ps(j)))
        end if
      end do
    end do
    end associate

  end function spin_spin_correlation


  pure function contract_2rdm(spinfree_2rdm, nActEl) result(spinfree_1rdm)

    !! Calculate the spinfree-1-RDM by tracing out one particle of the
    !! 3-index spinfree TwoRDM D(x1, x1', x2, x2). For debug purposes only.

    real(dp), intent(in) :: spinfree_2rdm(:,:,:,:)
    integer, intent(in) :: nelec

    real(dp), allocatable :: spinfree_1rdm(:,:)
    integer :: x, p, q

    mma_allocate(spinfree_1rdm(nActEl, nActEl))
    spinfree_1rdm(:,:) = 0.0_dp

    do x = 1, size(spinfree_2rdm, dim=1)
      do q = 1, size(spinfree_2rdm, dim=1)
        do p = 1, size(spinfree_2rdm, dim=1)
          spinfree_1rdm(p,q) = spinfree_1rdm(p,q) + spinfree_2rdm(p,q,x,x)
        end do
      end do
    end do

    spinfree_1rdm(:,:) = spinfree_1rdm(:,:) * 1.0_dp / (nelec - 1)

  end function contract_2rdm


  pure function decompress_symmetrized_2rdm(psmat, pamat) &
    result(spinfree_2rdm)
    ! tested in Python

    !! Convert a pair of PSMAT and PAMAT from linearised Molcas format
    !! into a full four-index spinfree-TwoRDM. For debug purposes only.

    real(dp), intent(in) :: pamat(:), psmat(:)
      !! (anti)symmetrised TwoRDM in the Molcas format.

    real(dp), allocatable :: spinfree_2rdm(:,:,:,:)

    mma_allocate(spinfree_2rdm(norbs, norbs, norbs, norbs))
    spinfree_2rdm(:,:,:,:) = 0.0_dp

    do pqrs = 1, size(psmat)
      call two_el_idx(pqrs, p, q, r, s)
      if (r == s) n_rs = 1
      if (r /= s) n_rs = 2
      twordm(p,q,r,s) = 2.0_dp/n_rs * (psmat(pqrs) + pamat(pqrs))
      twordm(r,s,p,q) = twordm(p,q,r,s)
      twordm(q,p,s,r) = 2.0_dp/n_rs * (psmat(pqrs) + pamat(pqrs))
      twordm(s,r,q,p) = twordm(q,p,s,r)
      twordm(q,p,r,s) = 2.0_dp/n_rs * (psmat(pqrs) - pamat(pqrs))
      twordm(r,s,q,p) = twordm(q,p,r,s)
      twordm(p,q,s,r) = 2.0_dp/n_rs * (psmat(pqrs) - pamat(pqrs))
      twordm(s,r,p,q) = twordm(p,q,s,r)
    end do

  end function decompress_symmetrized_2rdm

end module spin_correlation
