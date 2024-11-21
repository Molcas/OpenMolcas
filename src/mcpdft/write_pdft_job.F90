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
! 2023, Matthew R. Hennefarth - modified to modern fortran             *
!***********************************************************************

module write_pdft_job
  implicit none
  private

  public :: writejob

contains
  !> @brief Write energy and rotation matrix to JobIph or h5
  !>
  !> @author Matthew R. Hennefarth
  !>
  !> @param[in] e_states final PDFT energy for each state
  !> @param[in] nstates number of states
  !> @param[in] si_pdft Optional orthonormal eigenvectors to diagonalize effective Hamiltonian
  subroutine writejob(e_states,nstates,si_pdft)
    use definitions,only:iwp,wp
    use constants,only:zero
    use general_data,only:mxroot,mxiter

    integer(kind=iwp),intent(in) :: nstates
    real(kind=wp),intent(in) :: e_states(nstates)
    real(kind=wp),optional,intent(in) :: si_pdft(nstates,nstates)

    real(kind=wp) :: energy(mxroot*mxiter)
    integer(kind=iwp) :: i,j ! Dummy index variables for loops

    ! get energies
    energy = zero
    Do i = 1,mxIter
      Do j = 1,nstates
        energy(mxRoot*(i-1)+j) = e_states(j)
      Enddo
    Enddo

    call save_energies(energy)

    if(present(si_pdft)) then
      call save_ci(si_pdft,nstates)
    endif
  endsubroutine writejob

  !> @brief Save energies to wavefunction file
  !>
  !> @author Matthew R. Hennefarth
  !>
  !> @param[in] e_states state energies
  subroutine save_energies(e_states)
    use definitions,only:iwp,wp
    use general_data,only:jobiph
    use mcpdft_input,only:mcpdft_options
#ifdef _HDF5_
    use mh5,only:mh5_open_file_rw,mh5_open_dset,mh5_put_dset,mh5_close_file
#endif

#include "intent.fh"

    real(kind=wp),intent(_IN_) :: e_states(:)

    integer(kind=iwp) :: disk,adr19(15)
#ifdef _HDF5_
    integer(kind=iwp) :: refwfn_id,wfn_energy
#endif

    if(.not. mcpdft_options%is_hdf5_wfn) then
      adr19(:) = 0
      disk = 0
      call iDaFile(JOBIPH,2,adr19,15,disk)
      disk = adr19(6)
      call DDaFile(jobiph,1,e_states,size(e_states),disk)
#ifdef _HDF5_
    else
      refwfn_id = mh5_open_file_rw(mcpdft_options%wfn_file)
      wfn_energy = mh5_open_dset(refwfn_id,'ROOT_ENERGIES')
      call mh5_put_dset(wfn_energy,e_states)
      call mh5_close_file(refwfn_id)
#endif
    endif

  endsubroutine save_energies

  !> @brief Save MSPDFT final eigenvectors to wavefunction file
  !>
  !> @author Matthew R. Hennefarth
  !>
  !> @param[in] si_pdft Rotation matrix to final eigenstate basis
  !> @param[in] nstates number of states
  subroutine save_ci(si_pdft,nstates)
    use constants,only:zero,one
    use definitions,only:iwp,wp
    use general_data,only:jobiph
    use stdalloc,only:mma_allocate,mma_deallocate
    use mcpdft_input,only:mcpdft_options
#ifdef _HDF5_
    use mh5,only:mh5_open_file_rw,mh5_open_dset,mh5_put_dset, &
                  mh5_close_file,mh5_fetch_attr,mh5_fetch_dset
#endif

    integer(kind=iwp),intent(in) :: nstates
    real(kind=wp),intent(in) :: si_pdft(nstates,nstates)

    integer(kind=iwp) :: disk,ncon,state,ad19,dum(1),adr19(15)
    real(kind=wp),allocatable :: ci_rot(:,:),tCI(:,:)

#ifdef _HDF5_
    integer(kind=iwp) :: refwfn_id,wfn_cicoef
#endif

    ncon = 0

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

    call mma_allocate(tCI,nstates,ncon,label="tCI")
    call mma_allocate(ci_rot,nstates,ncon,label='CI Rot')

    if(.not. mcpdft_options%is_hdf5_wfn) then
      adr19(:) = 0
      ad19 = 0
      call iDaFile(JOBIPH,2,adr19,15,ad19)
      disk = adr19(4)
    endif

    do state = 1,nstates
      if(.not. mcpdft_options%is_hdf5_wfn) then
        call DDafile(jobiph,2,tCI(:,state),ncon,disk)
#ifdef _HDF5_
      else
        call mh5_fetch_dset(refwfn_id,'CI_VECTORS',tCI(:,state),[ncon,1],[0,state-1])
#endif
      endif
    enddo

    call dgemm_('n','n',nstates,ncon,nstates,one,si_pdft,nstates,tCI,ncon,zero,ci_rot,nstates)

    if(.not. mcpdft_options%is_hdf5_wfn) then
      disk = adr19(4)
      do state = 1,nstates
        call DDafile(jobiph,1,ci_rot(:,state),ncon,disk)
      enddo
#ifdef _HDF5_
    else
      wfn_cicoef = mh5_open_dset(refwfn_id,'CI_VECTORS')
      do state = 1,nstates
        call mh5_put_dset(wfn_cicoef,ci_rot(:,state),[ncon,1],[0,state-1])
      enddo
      call mh5_close_file(refwfn_id)
#endif
    endif

    call mma_deallocate(ci_rot)
    call mma_deallocate(tCI)

  endsubroutine save_ci

endmodule write_pdft_job
