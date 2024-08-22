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
!***********************************************************************
!                                                                      *
! 2023, Matthew R. Hennefarth - modified to modern fortran             *
!***********************************************************************

module write_pdft_job
  implicit none
  private

  public :: writejob

contains
  subroutine writejob(adr19,e_mspdft,si_pdft)
    ! Writes energy and rotation matrix (to final states) to
    ! either the jobiph or the h5 file.
    !
    ! Args:
    !   adr19: integer ndarray of shape (15)
    !       Holds info on location of where data is stored in jobiph
    !
    !   e_mspdft: ndarray of length lroots (optional)
    !       Array containin final MS-PDFT energies.
    !       Expected to be of length lroots (defined in rasscf.fh)
    !
    !   si_pdft: ndarray of length lroots*lroots (optional)
    !       Orthonormal eigenvectors of MS-PDFT in the intermediate
    !       state basis. Expected to be of length lroots*lroots.

    use definitions,only:wp
    implicit none

#include "WrkSpc.fh"
#include "rasdim.fh"
#include "rasscf.fh"
#include "general.fh"

    integer,intent(in) :: adr19(15)
    real(kind=wp),dimension(lroots),optional,intent(in) :: e_mspdft
    real(kind=wp),dimension(lroots**2),optional,intent(in) :: si_pdft
    real(kind=wp),dimension(mxroot*mxiter) :: energy
    real(kind=wp),dimension(lroots,lroots) :: U

    integer :: i,j ! Dummy index variables for loops

    ! get energies
    call dcopy_(mxRoot*mxIter,[0.0d0],0,energy,1)
    if(present(e_mspdft)) then
      Do i = 1,mxIter
        Do j = 1,lroots
          energy(mxRoot*(i-1)+j) = e_mspdft(j)
        Enddo
      Enddo
    else
      do i = 1,mxIter
        do j = 1,lroots
          energy(mxroot*(i-1)+j) = ener(j,1)
        enddo
      enddo
    endif

    call save_energies(adr19,energy)

    if(present(si_pdft)) then
      ! Move the rotated matrix into U variable
      do i = 1,lroots
        do j = 1,lroots
          U(j,i) = si_pdft(lroots*(i-1)+j)
        enddo
      enddo
      call save_ci(adr19,u)
    endif
  endsubroutine writejob

  subroutine save_energies(adr19,energy)
    ! Save the energies in the appropriate place (jobIPH or .h5 file)
    ! Args:
    !   adr19:
    !
    !   energy: ndarray of len mxroot*mxiter
    !     Final PDFT energies with zeros in the rest of the array.

    use definitions,only:wp
#ifdef _HDF5_
    use mh5,only:mh5_open_file_rw,mh5_open_dset,mh5_put_dset,mh5_close_file
#endif
    use mcpdft_input,only:mcpdft_options
    implicit none

    !for general.fh (needs mxSym)
#include "rasdim.fh"
    ! for jobiph
#include "general.fh"

    integer,dimension(15),intent(in) :: adr19
    integer :: disk
    real(kind=wp),dimension(:),intent(inout) :: energy
#ifdef _HDF5_
    integer :: refwfn_id,wfn_energy
#endif

    if(.not. mcpdft_options%is_hdf5_wfn) then
      disk = adr19(6)
      call DDaFile(jobiph,1,energy,size(energy),disk)
#ifdef _HDF5_
    else
      refwfn_id = mh5_open_file_rw(mcpdft_options%wfn_file)
      wfn_energy = mh5_open_dset(refwfn_id,'ROOT_ENERGIES')
      call mh5_put_dset(wfn_energy,energy(1))
      call mh5_close_file(refwfn_id)
#endif
    endif

  endsubroutine save_energies

  subroutine save_ci(adr19,U)
    ! Save the MS-PDFT final eigenvectors to either the jobIPH or .h5
    ! file.
    !
    ! Args:
    !   adr19:
    !
    !   U: ndarray of shape (lroots, lroots)
    !     Rotation matrix from intermediate state basis to final
    !     MS-PDFT eigenstate basis

    use definitions,only:wp
    use stdalloc,only:mma_allocate,mma_deallocate
#ifdef _HDF5_
    use mh5,only:mh5_open_file_rw,mh5_open_dset,mh5_put_dset, &
                  mh5_close_file,mh5_fetch_attr,mh5_fetch_dset
#endif
    use mcpdft_input,only:mcpdft_options
    implicit none

    ! for rasscf.fh ...
#include "rasdim.fh"
    ! for jobiph
#include "general.fh"
    integer,dimension(15),intent(in) :: adr19
    real(kind=wp),dimension(:,:),intent(in) :: U

    integer :: disk,ncon = 0,i,j,k
    integer,dimension(1) :: dum
    real(kind=wp),allocatable :: ci_rot(:),tCI(:)
    integer :: roots

#ifdef _HDF5_
    integer :: refwfn_id,wfn_cicoef
#endif

    roots = size(U,dim=1)

    if(.not. mcpdft_options%is_hdf5_wfn) then
      disk = 284 ! where does this number come from?
      call iDafile(jobiph,2,dum,1,disk)
      ncon = dum(1)

#ifdef _HDF5_
    else
      refwfn_id = mh5_open_file_rw(mcpdft_options%wfn_file)
      call mh5_fetch_attr(refwfn_id,'NCONF',ncon)
#endif
    endif

    call mma_allocate(ci_rot,ncon*roots,label='CI Rot')
    ci_rot = 0.0d0

    call mma_allocate(tCI,ncon,label="tCI")

    if(.not. mcpdft_options%is_hdf5_wfn) then
      disk = adr19(4)
    endif

    do i = 1,roots
      if(.not. mcpdft_options%is_hdf5_wfn) then
        call DDafile(jobiph,2,tCI,ncon,disk)
#ifdef _HDF5_
      else
        call mh5_fetch_dset(refwfn_id,'CI_VECTORS',tCI,[ncon,1],[0,i-1])
#endif
      endif

      do j = 1,roots
        do k = 1,ncon
          ci_rot((j-1)*ncon+k) = ci_rot((j-1)*ncon+k)+tCI(k)*U(i,j)
        enddo
      enddo
    enddo

    if(.not. mcpdft_options%is_hdf5_wfn) then
      disk = adr19(4)
      do i = 1,roots
        call DDafile(jobiph,1,ci_rot((i-1)*ncon+1),ncon,disk)
      enddo
#ifdef _HDF5_
    else
      wfn_cicoef = mh5_open_dset(refwfn_id,'CI_VECTORS')
      do i = 1,roots
        call mh5_put_dset(wfn_cicoef,ci_rot(ncon*(i-1)+1:ncon*i),[ncon,1],[0,i-1])
      enddo
      call mh5_close_file(refwfn_id)
#endif
    endif

    call mma_deallocate(ci_rot)
    call mma_deallocate(tCI)

  endsubroutine save_ci

endmodule write_pdft_job
