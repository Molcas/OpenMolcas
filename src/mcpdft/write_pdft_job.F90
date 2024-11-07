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
  subroutine writejob(e_pdft,si_pdft)
    ! Writes energy and rotation matrix (to final states) to
    ! either the jobiph or the h5 file.
    !
    ! Args:
    !   e_pdft: ndarray of length lroots (optional)
    !       Array containin final MS-PDFT energies.
    !       Expected to be of length lroots (defined in
    !       rasscf_global.F90)
    !
    !   si_pdft: ndarray of length lroots*lroots (optional)
    !       Orthonormal eigenvectors of MS-PDFT in the intermediate
    !       state basis. Expected to be of length lroots*lroots.

    use definitions,only:wp
    use constants,only:zero
    use rasscf_global,only:lRoots

    implicit none

#include "rasdim.fh"
#include "general.fh"

    real(kind=wp),dimension(lroots),intent(in) :: e_pdft
    real(kind=wp),dimension(lroots**2),optional,intent(in) :: si_pdft
    real(kind=wp),dimension(mxroot*mxiter) :: energy
    real(kind=wp),dimension(lroots,lroots) :: U

    integer :: i,j ! Dummy index variables for loops

    ! get energies
    energy = zero
    Do i = 1,mxIter
      Do j = 1,lroots
        energy(mxRoot*(i-1)+j) = e_pdft(j)
      Enddo
    Enddo

    call save_energies(energy)

    if(present(si_pdft)) then
      ! Move the rotated matrix into U variable
      do i = 1,lroots
        do j = 1,lroots
          U(j,i) = si_pdft(lroots*(i-1)+j)
        enddo
      enddo
      call save_ci(u)
    endif
  endsubroutine writejob

  subroutine save_energies(energy)
    ! Save the energies in the appropriate place (jobIPH or .h5 file)
    ! Args:
    !   energy: ndarray of len mxroot*mxiter
    !     Final PDFT energies with zeros in the rest of the array.

    use definitions,only:wp
#ifdef _HDF5_
    use mh5,only:mh5_open_file_rw,mh5_open_dset,mh5_put_dset,mh5_close_file
#endif
    use mcpdft_input,only:mcpdft_options
    implicit none

    real(kind=wp),dimension(:),intent(inout) :: energy

    !for general.fh (needs mxSym)
#include "rasdim.fh"
    ! for jobiph
#include "general.fh"

    integer,dimension(15) :: adr19
    integer :: disk,ad19
#ifdef _HDF5_
    integer :: refwfn_id,wfn_energy
#endif

    if(.not. mcpdft_options%is_hdf5_wfn) then
      adr19(:) = 0
      ad19 = 0
      call iDaFile(JOBIPH,2,adr19,15,ad19)
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

  subroutine save_ci(U)
    ! Save the MS-PDFT final eigenvectors to either the jobIPH or .h5
    ! file.
    !
    ! Args:
    !   U: ndarray of shape (lroots, lroots)
    !     Rotation matrix from intermediate state basis to final
    !     MS-PDFT eigenstate basis

    use constants,only:zero,one
    use definitions,only:iwp,wp
    use stdalloc,only:mma_allocate,mma_deallocate
#ifdef _HDF5_
    use mh5,only:mh5_open_file_rw,mh5_open_dset,mh5_put_dset, &
                  mh5_close_file,mh5_fetch_attr,mh5_fetch_dset
#endif
    use mcpdft_input,only:mcpdft_options
    implicit none

#include "rasdim.fh"
    ! for jobiph
#include "general.fh"
    integer(kind=iwp),dimension(15) :: adr19
    real(kind=wp),dimension(:,:),intent(in) :: U

    integer(kind=iwp) :: disk,ncon = 0,i,ad19
    integer(kind=iwp),dimension(1) :: dum
    real(kind=wp),allocatable :: ci_rot(:,:),tCI(:,:)
    integer(kind=iwp) :: roots

#ifdef _HDF5_
    integer(kind=iwp) :: refwfn_id,wfn_cicoef
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

    call mma_allocate(tCI,roots,ncon,label="tCI")
    call mma_allocate(ci_rot,roots,ncon,label='CI Rot')

    if(.not. mcpdft_options%is_hdf5_wfn) then
      adr19(:) = 0
      ad19 = 0
      call iDaFile(JOBIPH,2,adr19,15,ad19)
      disk = adr19(4)
    endif

    do i = 1,roots
      if(.not. mcpdft_options%is_hdf5_wfn) then
        call DDafile(jobiph,2,tCI(i,:),ncon,disk)
#ifdef _HDF5_
      else
        call mh5_fetch_dset(refwfn_id,'CI_VECTORS',tCI(i,:),[ncon,1],[0,i-1])
#endif
      endif
    enddo

    call dgemm_('n','n',roots,ncon,roots,one,U,roots,tCI,ncon,zero,ci_rot,roots)

    if(.not. mcpdft_options%is_hdf5_wfn) then
      disk = adr19(4)
      do i = 1,roots
        call DDafile(jobiph,1,ci_rot(i,:),ncon,disk)
      enddo
#ifdef _HDF5_
    else
      wfn_cicoef = mh5_open_dset(refwfn_id,'CI_VECTORS')
      do i = 1,roots
        call mh5_put_dset(wfn_cicoef,ci_rot(i,:),[ncon,1],[0,i-1])
      enddo
      call mh5_close_file(refwfn_id)
#endif
    endif

    call mma_deallocate(ci_rot)
    call mma_deallocate(tCI)

  endsubroutine save_ci

endmodule write_pdft_job
