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
* Copyright (C) 2020, Chen Zhou                                        *
************************************************************************
      subroutine writejobms(iadr19,L1,L2)
      use definitions, only: wp
      use stdalloc, only: mma_allocate, mma_deallocate
#ifdef _HDF5_
      use mh5, only: mh5_open_file_rw, mh5_open_dset,
     &               mh5_put_dset, mh5_close_file, mh5_fetch_attr,
     &               mh5_fetch_dset
      use write_pdft_job, only: save_ci,
     &                          save_energies
#else
      use write_pdft_job, only: save_ci, save_energies
#endif
      implicit none

#include "WrkSpc.fh"
#include "wadr.fh"
#include "rasdim.fh"
#include "rasscf.fh"
#include "general.fh"

      integer :: i,j,L1,L2
      integer :: iadr19(15)
      real(kind=wp), allocatable :: energy(:)

      ! get energies
      call mma_allocate(energy, mxroot*mxiter, label="Energy")
      call dcopy_(mxRoot*mxIter,[0.0d0],0,energy,1)
      Do i = 1,mxIter
        Do j = 1,lroots
          energy(mxRoot*(i-1)+j) = work(L1+j-1)
        End do
      End do

      call save_energies(iadr19, energy)
      call mma_deallocate(energy)

      call save_ci(iadr19, l2)

      End

      subroutine writejob(iadr19)
      use definitions, only: wp
      use stdalloc, only: mma_allocate, mma_deallocate
      use write_pdft_job, only: save_energies
      implicit none

    ! For mxRoot and mxIter
#include "rasdim.fh"
    ! For ener
#include "rasscf.fh"

      integer, dimension(15), intent(in) :: iadr19
      integer :: i, j
      real(kind=wp), allocatable :: energy(:)

      call mma_allocate(energy, mxRoot*mxIter,label="energy")
      call dcopy_(mxRoot*mxIter, [0.0d0], 0, energy, 1)
      do i=1, mxIter
        do j=1, lroots
          energy(mxroot*(i-1)+j) = ener(j,1)
        end do
      end do

      call save_energies(iadr19, energy)

      call mma_deallocate(energy)

      End
