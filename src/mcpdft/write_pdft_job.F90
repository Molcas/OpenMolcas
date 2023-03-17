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
! Copyright (C) 2020, Chen Zhou                                        *
! Copyright (C) 2023, Matthew R. Hennefarth                            *
!***********************************************************************

module write_pdft_job
  implicit none
  private
  integer :: iwjob

  ! read from / write to HDF5 file
  logical :: hasHDF5ref=.false.

  ! reference wave function is of MPS type (aka "DMRG wave function")
  logical :: hasMPSref=.false.

  public :: iwjob, hasHDF5ref, hasMPSref, save_energies, save_ci

  contains
  subroutine save_energies(adr19, energy)
    use definitions, only: wp
#ifdef _HDF5_
    use mh5, only: mh5_open_file_rw, mh5_open_dset, mh5_put_dset, mh5_close_file
#endif
    implicit none

    !for general.fh (needs mxSym)
#include "rasdim.fh"
    ! for jobiph
#include "general.fh"

    integer, dimension(15), intent(in) :: adr19
    integer :: disk
    real(kind=wp), dimension(:), intent(in) :: energy
#ifdef _HDF5_
    integer :: refwfn_id, wfn_energy
#endif

    if (.not.hasHDF5ref) then
      disk = adr19(6)
      call DDaFile(jobiph, 1, energy, size(energy), disk)
#ifdef _HDF5_
    else
      if (hasMPSref) then
        refwfn_id = mh5_open_file_rw(StartOrbFile)
      else
        refwfn_id = mh5_open_file_rw('RASWFN')
      end if
      wfn_energy = mh5_open_dset(refwfn_id, 'ROOT_ENERGIES')
      call mh5_put_dset(wfn_energy, energy(1))
      call mh5_close_file(refwfn_id)
#endif
    end if

  end subroutine save_energies

  subroutine save_ci(adr19, l2)
    use definitions, only: wp
    use stdalloc, only: mma_allocate, mma_deallocate
#ifdef _HDF5_
      use mh5, only: mh5_open_file_rw, mh5_open_dset, mh5_put_dset, &
                     mh5_close_file, mh5_fetch_attr,mh5_fetch_dset
#endif
    implicit none

    ! for work
#include "WrkSpc.fh"
    ! for rasscf.fh ...
#include "rasdim.fh"
    ! for lroots
#include "rasscf.fh"
    ! for jobiph
#include "general.fh"
    integer, dimension(15), intent(in) :: adr19
    integer, intent(in) :: l2

    integer :: disk, ncon, i, j, k
    integer, dimension(1) :: dum
    real(kind=wp), dimension(lroots, lroots) :: U
    real(kind=wp), allocatable :: ci_rot(:), tCI(:)

#ifdef _HDF5_
    integer :: refwfn_id, wfn_cicoef
#endif

    if (.not.hasHDF5ref) then
      disk = 284 ! where does this number come from?
      call iDafile(jobiph, 2, dum, 1, disk)
      ncon = dum(1)

#ifdef _HDF5_
    else
      if (hasMPSref) then
        refwfn_id = mh5_open_file_rw(StartOrbFile)
      else
        refwfn_id = mh5_open_file_rw('RASWFN')
      end if
      call mh5_fetch_attr(refwfn_id, 'NCONF', ncon)
#endif
    end if

    do i=1, lroots
      do j=1, lroots
        U(j,i) = work(L2 + (i-1)*lroots+j-1)
      end do
    end do

    call mma_allocate(ci_rot, ncon*lroots, label='CI Rot')
    call dcopy_(ncon*lroots, [0.0d0], 0, ci_rot, 1)

    call mma_allocate(tCI, ncon, label="tCI")

    if (.not.hasHDF5ref) then
      disk = adr19(4)
    end if

    do i=1, lroots
      if(.not.hasHDF5ref) then
        call DDafile(jobiph, 2, tCI, ncon, disk)
#ifdef _HDF5_
      else
        call mh5_fetch_dset(refwfn_id, 'CI_VECTORS', tCI, [ncon,1], [0,i-1])
#endif
      end if

      do j=1, lroots
        do k=1, ncon
          ci_rot((j-1)*ncon + k) = ci_rot((j-1)*ncon + k) + tCI(k)*U(i,j)
        end do
      end do
    end do

    if (.not.hasHDF5ref) then
      disk = adr19(4)
      do i=1, lroots
        call DDafile(jobiph, 1, ci_rot((i-1)*ncon+1), ncon, disk)
      end do
#ifdef _HDF5_
    else
      wfn_cicoef = mh5_open_dset(refwfn_id, 'CI_VECTORS')
      do i=1, lroots
        call mh5_put_dset(wfn_cicoef, ci_rot(ncon*(i-1)+1:ncon*i), [ncon, 1], [0,i-1])
      end do
      call mh5_close_file(refwfn_id)
#endif
    end if

    call mma_deallocate(ci_rot)
    call mma_deallocate(tCI)

  end subroutine save_ci

end module write_pdft_job