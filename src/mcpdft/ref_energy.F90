!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!***********************************************************************

subroutine ref_energy(mcscf_energy,nstates)
  use definitions,only:iwp,wp,u6
  use constants,only:zero
  use stdalloc,only:mma_allocate,mma_deallocate
  use printlevel,only:usual
  use mcpdft_output,only:iprglb
  use mcpdft_input,only:mcpdft_options
  use mspdft,only:heff
  use general_data,only:jobiph,mxiter,mxroot
#ifdef _HDF5_
  use mh5,only:mh5_open_file_r,mh5_close_file,mh5_fetch_dset,mh5_exists_dset
#endif

  implicit none

  integer(kind=iwp),intent(in) :: nstates
  real(kind=wp),intent(out) :: mcscf_energy(nstates)

  integer(kind=iwp) :: state,iad19,disk,nmaybe,i,it,iadr19(15)
  real(kind=wp) :: e,aemax
  real(kind=wp),allocatable :: elist(:,:)

#ifdef _HDF5_
  integer(kind=iwp) :: refwfn_id
#endif

  if(mcpdft_options%mspdft) then
    if(iprglb >= usual) then
      write(u6,*) 'Reference MC-SCF energies taken from diagonal elements of'
      write(u6,*) 'effective Hamiltonian'
    endif
    do state = 1,nstates
      mcscf_energy(state) = heff(state,state)
    enddo
  else
    if(iprglb >= usual) then
      write(u6,*) 'Reference MC-SCF energies taken from ',mcpdft_options%wfn_file
    endif
    if(.not. mcpdft_options%is_hdf5_wfn) then
      iadr19(:) = 0
      iad19 = 0
      call iDaFile(JOBIPH,2,iadr19,15,iad19)
      disk = iadr19(6)

      call mma_allocate(elist,mxroot,mxiter,label='EList')
      elist = zero
      call DDaFile(JOBIPH,2,elist,mxroot*mxiter,disk)

      nmaybe = 0
      do it = 1,mxiter
        aemax = zero
        do i = 1,mxroot
          e = elist(i,it)
          aemax = max(aemax,abs(e))
        enddo
        if(abs(aemax) <= 1.0D-12) then
          exit
        endif
        nmaybe = it
      enddo

      mcscf_energy = elist(:,nmaybe)
      do state = 1,nstates
        mcscf_energy(state) = elist(state,nmaybe)
      enddo

      call mma_deallocate(elist)

#ifdef _HDF5_
    else
      refwfn_id = mh5_open_file_r(mcpdft_options%wfn_file)
      if(.not. mh5_exists_dset(refwfn_id,'ROOT_ENERGIES')) then
        write(u6,*) 'The HDF5 ref file does not contain ROOT_ENERGIES'
        write(u6,*) 'Fatal error, the calculation will stop now.'
        call Abend()
      endif
      write(u6,*) 'Loading mcscf energy from hdf5 file!'
      call mh5_fetch_dset(refwfn_id,'ROOT_ENERGIES',mcscf_energy)
      call mh5_close_file(refwfn_id)
#endif
    endif
  endif

endsubroutine
