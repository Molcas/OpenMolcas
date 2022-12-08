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

!> Read information provided on runfile/jobiph files and
!> initialize pt2wfn file
subroutine pt2init(refwfn_in)

#ifdef _DMRG_
use qcmaquis_info, only: qcm_group_names
use qcmaquis_interface, only: qcmaquis_interface_measure_and_save_trans3rdm, qcmaquis_interface_get_3rdm_elements, &
                              qcmaquis_interface_get_4rdm_elements
use qcmaquis_interface_utility_routines, only: str
#endif
#ifdef _HDF5_QCM_
use hdf5_utils, only: hdf5_init, hdf5_close, hdf5_exit, hdf5_get_data, hdf5_open
use mh5, only: mh5_is_hdf5
#endif
use refwfn, only: refwfn_init, refwfn_info, refwfn_data, refwfn_close
use nevpt2_cfg, only: curr_dir, do_cholesky, igelo, molcas_project, MultGroup, nr_active_electrons, nr_frozen_orb, nr_states, nspin
use info_state_energy, only: e, init_energies                                                     ! energies
use info_orbital_space, only: datadim, file_id, ijklname, inforb_molcas, initialize_inforb_molcas ! orbital specifications read from JobIph
use nevpt2wfn, only: nevpt2wfn_init, nevpt2wfn_data
use stdalloc, only: mma_allocate, mma_deallocate
use Definitions, only: wp, iwp, u6
use MxDM

implicit none
character(len=*), intent(in) :: refwfn_in
character(len=:), allocatable :: refwfnfile
integer(kind=iwp) :: istate, ii, j, nDiff, nishprev, nfroprev
integer(kind=iwp), allocatable :: nCore_local(:)
real(kind=wp), allocatable :: readbuf(:,:)
#include "caspt2.fh"

! Save current directory into the CurrDir string
call GetEnvF('CurrDir',curr_dir)
call GetEnvF('Project',molcas_project)

! call the Molcas routine to check whether we're using Cholesky
call DecideOnCholesky(do_cholesky)
! read some stuff from RunFile, just as CASPT2 does
call Get_iScalar('nSym',nSym)
! Symmetry with Cholesky decomposition is not supported yet!
if ((nSym > 1) .and. do_cholesky) then
  call WarningMessage(1,'Symmetry with Cholesky decomposition is not supported yet!')
  call Quit_OnUserError()
end if
write(u6,*) 'Cholesky Decomposition: ',merge('Enabled ','Disabled',do_cholesky)

!> set up AO basis information (needed for post-NEVPT2 information transfer)
call Get_iArray('nBas',nBas,nSym)
nbast = sum(nbas(1:nsym))
nbsqt = sum(nbas(1:nsym)**2)

!> TODO: Read more stuff, such as the information about the active
!> space from the JobIph file, or at least add some consistency checks

!> open files in analogy to CASPT2
LUONEM = 16
call DANAME_wa(LUONEM,'MOLONE')

! We can only open HDF5 files, so if no HDF5 file is specified in the input
! or JOBIPH is not a HDF5 file, try to find the corresponding HDF5 file,
! if not found, exit

#ifdef _HDF5_

#ifdef _WARNING_WORKAROUND_
! fix a compiler warning/error about possibly uninitialized variables
allocate(character(len=0) :: refwfnfile)
#endif
refwfnfile = trim(refwfn_in)
if (.not. mh5_is_hdf5(refwfnfile)) then
  ! try $Project.dmrgscf.h5
  refwfnfile = trim(molcas_project)//'.dmrgscf.h5'
  if (.not. mh5_is_hdf5(refwfnfile)) then
    ! try $Project.rasscf.h5
    refwfnfile = trim(molcas_project)//'.rasscf.h5'
    if (.not. mh5_is_hdf5(refwfnfile)) then
      call WarningMessage(1,'Cannot find a HDF5 file with the reference wavefunction. Make sure that file '// &
                          trim(molcas_project)//'.rasscf.h5 or .dmrgscf.h5 exists')
      call Quit_OnUserError()
    end if
  end if
end if
#else
call WarningMessage(1,'Please compile OpenMolcas with HDF5 support for NEVPT2 to work')
call Quit_OnUserError()
#endif
call refwfn_init(refwfnfile)
call refwfn_info()
call refwfn_data()
call refwfn_close()

!> fill nevpt2 confuguration variables from caspt2.fh commons
!> ----------------------------------------------------------

!> check if nr_states has been requested as 'all'
if (nr_states == 0) then
  ! using standard allocate and deallocate because MultGroup%State
  ! is deallocated somewhere in the external library
  if (allocated(MultGroup%State)) deallocate(MultGroup%State)
  !> nstate from common block in caspt2.fh
  nr_states = nstate
  allocate(MultGroup%State(nr_states))
  do istate=1,nr_states
    MultGroup%State(istate) = istate
  end do
end if

!> spin & number of active electrons
nspin = ispin
nr_active_electrons = nactel

write(u6,'(/a)') ' Wavefunction parameters for NEVPT2'
write(u6,'(a )') ' ----------------------------------'
write(u6,'(a,i4)') ' Number of active electrons ....... ',nr_active_electrons
write(u6,'(a,i4)') ' Spin ............................. ',nspin

! read the info about orbital spaces and initialise
! the NEVPT2 inforb array

! Read orbital specifications and store them in the inforb_molcas
! variable from info_orbital_space

! nish, nash and nssh are in caspt2.fh (Common /INPI/ and have been
! read by the refwfn module
call initialize_inforb_molcas(nSym)

if (maxval(nfro) /= 0) then
  ! Leon: a workaround if frozen orbitals in MCSCF have been detected:
  ! MCSCF frozen orbitals MUST count as active here and then added
  ! to igelo(:) later, otherwise the orbital numbers do not match!
  write(u6,*) 'Frozen orbitals in MCSCF detected.'
  inforb_molcas%nish(1:nSym) = nfro(1:nSym)+nish(1:nSym)
else
  inforb_molcas%nish(1:nSym) = nish(1:nSym)
end if

! Read frozen orbitals from ijkl.h5 files produced by MOTRA
! we need this information now to calculate the number of frozen
! orbitals correctly

! warning! are we sure that nsym here is the same as nsym in hdf5?

! is norb from MOTRA the same as nBas or not? or is it nBas-nfro-ndel?
! In any case, we're reading the value from MOTRA just to be consistent
#ifdef _HDF5_QCM_
call hdf5_init()
call hdf5_open(ijklname,file_id(2)) ! open ijkl.h5
datadim(1) = nsym
call mma_allocate(readbuf,nsym,3,label='readbuf')
readbuf(:,:) = -1
call hdf5_get_data(file_id(2),'norb  ',datadim,readbuf(1,1))
call hdf5_get_data(file_id(2),'nfro  ',datadim,readbuf(1,2))
call hdf5_get_data(file_id(2),'ndel  ',datadim,readbuf(1,3))
inforb_molcas%norb(1:nSym) = nint(readbuf(1:nSym,1))
inforb_molcas%nfro(1:nSym) = nint(readbuf(1:nSym,2))
inforb_molcas%ndel(1:nSym) = nint(readbuf(1:nSym,3))
call mma_deallocate(readbuf)

call hdf5_close(file_id(2)) ! close ijkl.h5
call hdf5_exit()
#else
! Should never be the case!
call WarningMessage(1,'HDF5 QCMaquis interface not enabled, cannot continue!')
call Abend()
#endif
! Detect frozen orbitals from reference wavefunction, if there was
! no input wrt frozen orbitals

if (.not. allocated(igelo)) then
  if (nr_frozen_orb /= -1) then
    ! Autodetect frozen orbitals:

    ! If nr_frozen_orb is set to -1 in rdinput()
    ! it's to signal that # of frozen orbs has
    ! been forcibly set to 0
    ! no allocated frozen orbitals -> no frozen input
    ! hence try to guess frozen orbitals from reference wavefunction

    ! Correct the # of frozen orbitals just as it is done in CASPT2
    call mma_allocate(nCore_local,nSym,label='nCore_local')
    call Get_iArray('Non valence orbitals',nCore_local,nSym)
    do ii=1,nSym
      if (nCore_local(ii) > nFro(ii)) then
        nDiff = nCore_local(ii)-nFro(ii)
        nDiff = min(nDiff,nISh(ii))
        nFro(ii) = nFro(ii)+nDiff
        nISh(ii) = nISh(ii)-nDiff
      end if
    end do
    call mma_deallocate(nCore_local)

    ! Check if orbitals have been frozen in MOTRA
    if (maxval(inforb_molcas%nfro(1:nSym)) > 0) then
      do ii=1,nSym
        nFro(ii) = nFro(ii)-inforb_molcas%nfro(ii)
        if (nFro(ii) < 0) then
          call WarningMessage(1,'Warning: Additional frozen orbitals in MOTRA.')
          nFro(ii) = 0
        end if
      end do
    end if

    ! calculate total # of frozen orbitals
    nr_frozen_orb = 0
    do ii=1,nSym
      nr_frozen_orb = nr_frozen_orb+nFro(ii)
    end do

    ! using standard allocate because igelo
    ! is deallocated somewhere in the external library
    allocate(igelo(nr_frozen_orb))
    ! fill igelo by symmetry
    !! TODO: Check here for what has been frozen in MOTRA and subtract it (!!)
    do ii=1,nSym
      if (ii == 1) then
        nishprev = 0
        nfroprev = 0
      else
        nishprev = nfro(ii-1)+nish(ii-1)
        nfroprev = nfro(ii-1)
      end if
      do j=1,nfro(ii)
        igelo(nfroprev+j) = nishprev+j
      end do
    end do
  else
    ! Reset nr_frozen_orb back to 0, as we have been signalled with
    ! nr_frozen_orb = -1 from rdinput that frozen orbs have been forcibly
    ! reset to zero
    nr_frozen_orb = 0
  end if
else
  if (maxval(inforb_molcas%nfro(1:nSym)) > 0) then
    call WarningMessage(2,'Warning: Orbitals have been frozen in MOTRA and manual freeze specification has been provided in '// &
                        'NEVPT2 input. Make sure you are not double-counting the frozen orbitals!')
  end if
end if

inforb_molcas%nash(1:nSym) = nash(1:nSym)
inforb_molcas%nssh(1:nSym) = nssh(1:nSym)
inforb_molcas%nbas(1:nSym) = nbas(1:nSym)
inforb_molcas%nbast = nbast
inforb_molcas%ncmo = ncmo
inforb_molcas%nbsqt = nbsqt

if (maxval(inforb_molcas%nfro(1:nSym)) > 0) then
  write(u6,'(a,8(18i4))') ' Frozen orbitals from MOTRA ....... ',(inforb_molcas%nfro(ii),ii=1,nSym)
end if

write(u6,'(a,8(18i4))') ' Inactive orbitals ................ ',(inforb_molcas%nish(ii),ii=1,nSym)
write(u6,'(a,8(18i4))') ' Active orbitals .................. ',(inforb_molcas%nash(ii),ii=1,nSym)
write(u6,'(a,8(18i4))') ' Secondary orbitals ............... ',(inforb_molcas%nssh(ii),ii=1,nSym)

if (maxval(inforb_molcas%ndel(1:nSym)) > 0) then
  write(u6,'(a,8(18i4))') ' Deleted orbitals from MOTRA..... ',(inforb_molcas%ndel(ii),ii=1,nSym)
end if

write(u6,'(a/)') ' ----------------------------------'

call init_energies(nr_states)

write(u6,'(a)') ' Energies of zeroth-order DMRG wavefunction(s)'
write(u6,'(a)') ' ---------------------------------------------'
do istate=1,nr_states
  !> copy reference energies
  e(istate) = refene(MultGroup%State(istate))
  write(u6,'(a,i4,a,f18.8)') ' State ...',MultGroup%State(istate),' ... Energy = ',e(istate)
end do
write(u6,'(a/)') ' ---------------------------------------------'

#ifdef _DMRG_
if (allocated(qcm_group_names)) then
  write(u6,'(a)') ' DMRG wavefunction data will be read from'
  write(u6,'(a)') ' ----------------------------------------'
  if (.not. allocated(MultGroup%h5_file_name)) call mma_allocate(MultGroup%h5_file_name,nr_states,label='h5_file_name')
  MultGroup%h5_file_name = ''
  do istate=1,nr_states
    MultGroup%h5_file_name(istate) = trim(qcm_group_names(1)%states(MultGroup%State(istate)))
    !> copy reference wfn file names
    write(u6,'(a,i4,a,a)') ' State ...',MultGroup%State(istate),' .......................... ', &
                           trim(qcm_group_names(1)%states(MultGroup%State(istate)))
  end do
  write(u6,'(a/)') ' ----------------------------------------'
end if
#endif

!> Create the PT2 wavefunction file. The reference file should not be active,
!> as it might be the same file (in which case it is overwritten).
call nevpt2wfn_init(.true.)
call nevpt2wfn_data()

end subroutine pt2init
