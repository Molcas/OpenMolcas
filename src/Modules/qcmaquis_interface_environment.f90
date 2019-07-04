!!  dmrg-interface-utils: interface to the Maquis DMRG program for various
!!                        quantum-chemistry program packages.
!!  Copyright 2013-2018 Leon Freitag, Erik Hedegaard, Sebastian Keller,
!!                      Stefan Knecht, Yingjin Ma, Christopher Stein
!!                      and Markus Reiher
!!                      Laboratory for Physical Chemistry, ETH Zurich
!!
!!  dmrg-interface-utils is free software: you can redistribute it and/or modify
!!  it under the terms of the GNU Lesser General Public License as published by
!!  the Free Software Foundation, either version 3 of the License, or
!!  (at your option) any later version.
!!
!!  dmrg-interface-utils is distributed in the hope that it will be useful,
!!  but WITHOUT ANY WARRANTY; without even the implied warranty of
!!  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
!!  GNU Lesser General Public License for more details.
!!
!!  You should have received a copy of the GNU Lesser General Public License
!!  along with dmrg-interface-utils. If not, see <http://www.gnu.org/licenses/>.

module qcmaquis_interface_environment

! stefan: transfer external MOLCAS settings to DMRG environment variables

  use qcmaquis_interface_cfg
  use qcmaquis_interface_utility_routines

  implicit none


  public initialize_dmrg
  public finalize_dmrg
  public initialize_dmrg_rassi
  public print_dmrg_info
  public dump_dmrg_info
  public read_dmrg_info
  public set_dmrg_runtime_environment
  public save_dmrg_parameter_for_mclr

contains

  subroutine initialize_dmrg(                         &
     !>>>>>>>>>>>>>>>>>>>>>>>>>>>>   General Parameters  <<<<<<<<<<<<<<<<<<<<<<<<<<<!
                             nsym_molcas,             &
                             lsym_molcas,             &
                             nactel_molcas,           &
                             ispin_molcas,            &
                             nroots_molcas,           &
                             maxroot_molcas,          &
                             iroot_molcas,            &
                             nrs2_molcas,             &
                             LRras2_mclr_molcas,      &
                             tash_molcas,             &
                             weight_molcas,           &
                             thre_molcas,             &
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>   DMRG                 <<<<<<<<<<<<<<<<<<<<<<<<<<<!
                             dmrg_molcas,             &
                             init_DMRG,               &
                             initial_occ              &

    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>   MPI                  <<<<<<<<<<<<<<<<<<<<<<<<<<<!
#ifdef _MOLCAS_MPP_
                            ,nprocs_mpi,              &
                             myrank_mpi               &
#endif
                            )

    !> external MOLCAS variables
    integer, intent(in)                           :: nsym_molcas
    integer, intent(in)                           :: lsym_molcas
    integer, intent(in)                           :: nactel_molcas
    integer, intent(in)                           :: ispin_molcas
    integer, intent(in)                           :: nroots_molcas
    integer, intent(in)                           :: maxroot_molcas
    integer, intent(in)                           :: tash_molcas
    integer, intent(in), dimension(maxroot_molcas):: iroot_molcas
    integer, intent(in), dimension(nsym_molcas)   :: nrs2_molcas
    integer, intent(in), dimension(nsym_molcas)   :: LRras2_mclr_molcas
    integer, optional, intent(in)                 :: initial_occ(tash_molcas,nroots_molcas)
#ifdef _MOLCAS_MPP_
    integer, intent(in)                           :: nprocs_mpi
    integer, intent(in)                           :: myrank_mpi
#endif

    real*8 , intent(in)                           :: thre_molcas
    real*8 , intent(in), dimension(maxroot_molcas):: weight_molcas

    logical, intent(in)                           :: dmrg_molcas

    character*20, intent(in)                      :: init_DMRG

    !> local variables
    integer              :: ierr, i, lp, lcd
    integer              :: nasht
    integer              :: nproc_DMRG
    integer              :: multd2h(8,8)
    integer, allocatable :: nash(:)
    character(len=20)    :: check_nr_cores
    character(len=36)    :: check_fc_offset
    character(len=1695)  :: currdir
    character(len=300)   :: project

    !> control flags
         doDMRG = dmrg_molcas
    !> Threshold for energy
    E_threshold = thre_molcas


    if(.not. doDMRG) return

    multd2h(1,:) = (/1,2,3,4,5,6,7,8/)
    multd2h(2,:) = (/2,1,4,3,6,5,8,7/)
    multd2h(3,:) = (/3,4,1,2,7,8,5,6/)
    multd2h(4,:) = (/4,3,2,1,8,7,6,5/)
    multd2h(5,:) = (/5,6,7,8,1,2,3,4/)
    multd2h(6,:) = (/6,5,8,7,2,1,4,3/)
    multd2h(7,:) = (/7,8,5,6,3,4,1,2/)
    multd2h(8,:) = (/8,7,6,5,4,3,2,1/)

    dmrg_symmetry = type_symmetry(                 &
                                  nsym_molcas,     &
                                  multd2h(1:8,1:8) &
                                 )

    !> default # CPUS: try first QCMaquis_CPUS and then OMP_NUM_THREADS if the first one is not available.
    nproc_dmrg = 1
    call getenv("QCMaquis_CPUS",check_nr_cores)
    if(check_nr_cores /= "")then
      read(check_nr_cores,'(I5)') nproc_dmrg
    else
      call getenv("OMP_NUM_THREADS",check_nr_cores)
      if(check_nr_cores /= "") read(check_nr_cores,'(I5)') nproc_dmrg
    end if

    dmrg_setup    = type_setup             (                                &
                                            init_DMRG(1:2),                 &
                                            nproc_DMRG                      &
                                           )
    !> set prefix for QCMaquis results/checkpoint files
    call getenv("Project",project)
    call getenv("CurrDir",currdir)

    dmrg_file%prefix(1:1997) = ''
    lp = 0; lcd = 0
    lp                       = len_trim(project)
    lcd                      = len_trim(currdir)
    if(trim(project) /= '' .and. trim(currdir) /= '')then
      dmrg_file%prefix(1:lp+lcd+2) = trim(currdir)//'/'//trim(project)//'.'
    else
      dmrg_file%prefix(1:2) = './'
    end if

    !> set offset for QCMaquis result/checkpoint file counting (useful if we run several RASSCF jobs in a single shot)
    dmrg_file%offset = 0
    call getenv("QCMFCO",check_fc_offset)
    if(check_fc_offset /= "")then
      read(check_fc_offset,'(I5)') dmrg_file%offset
    end if

    allocate(dmrg_state%iroot(maxroot_molcas), stat=ierr); if( ierr /= 0 )stop ' Error in allocation: iroot(:)'
    allocate(dmrg_state%weight(maxroot_molcas), stat=ierr); if( ierr /= 0 )stop ' Error in allocation: weight(:)'

    dmrg_state    = type_state             (                                &
                                            lsym_molcas,                    &
                                            nactel_Molcas,                  &
                                            iSpin_molcas,                   &
                                            iSpin_molcas-1,                 &
                                            nroots_molcas,                  &
                                            maxroot_molcas,                 &
                                            iroot_Molcas(1:maxroot_molcas), &
                                            weight_Molcas(1:maxroot_molcas) &
                                           )

    allocate(nash(nsym_molcas))

    nasht = 0
    do i = 1, nsym_molcas
      nash(i)    = nrs2_molcas(i)
      nasht      = nasht + nash(i)
    end do

    ! Leon: Initialise dmrg_external%norb already here
    dmrg_external%norb   = nasht

    allocate(dmrg_orbital_space%initial_occ(nasht,nroots_molcas), stat=ierr); if( ierr /= 0 ) &
    stop ' Error in allocation: initial_occ(:,:)'
    dmrg_orbital_space = type_orbital_space(                                     &
                                            nash(1:nsym_molcas),                 &
                                            LRras2_mclr_molcas(1:nsym_molcas),   &
                                            initial_occ(1:nasht,1:nroots_molcas) &
                                           )
    deallocate(nash)

    allocate(dmrg_energy%dmrg_state_specific(maxroot_molcas), stat=ierr); if( ierr /= 0 ) &
    stop ' Error in allocation: dmrg_state_specific(:)'
    dmrg_energy%dmrg_state_specific = 0.0d0

    allocate(dmrg_energy%num_sweeps(maxroot_molcas), stat=ierr); if( ierr /= 0 ) &
    stop ' Error in allocation: num_sweeps(:)'
    dmrg_energy%num_sweeps     = 0

    allocate(dmrg_energy%num_sweeps_old(maxroot_molcas), stat=ierr); if( ierr /= 0 ) &
    stop ' Error in allocation: num_sweeps_old(:)'
    dmrg_energy%num_sweeps_old = 0

    allocate(dmrg_energy%max_truncW(maxroot_molcas), stat=ierr); if( ierr /= 0 ) &
    stop ' Error in allocation: max_truncW(:)'
    dmrg_energy%max_truncW     = 0
    allocate(dmrg_energy%max_truncW_old(maxroot_molcas), stat=ierr); if( ierr /= 0 ) &
    stop ' Error in allocation: max_truncW_old(:)'
    dmrg_energy%max_truncW_old     = 0

     allocate(dmrg_file%qcmaquis_parameter_file(maxroot_molcas), stat=ierr); if( ierr /= 0 ) &
     stop ' Error in allocation: qcmaquis_parameter_file(:)'
     dmrg_file%qcmaquis_parameter_file = ''

     allocate(dmrg_file%qcmaquis_checkpoint_file(maxroot_molcas), stat=ierr); if( ierr /= 0 ) &
     stop ' Error in allocation: qcmaquis_checkpoint_file(:)'
     dmrg_file%qcmaquis_checkpoint_file = ''

     allocate(dmrg_orbital_ordering%fiedler_order(maxroot_molcas), stat=ierr); if( ierr /= 0 ) &
     stop ' Error in allocation: fiedler_order(:)'
     dmrg_orbital_ordering%fiedler_order = ''

    !> initialize parallel settings from the host program
#ifdef _MOLCAS_MPP_
    dmrg_host_program_settings%nprocs = nprocs_mpi
    dmrg_host_program_settings%myrank = myrank_mpi
    if(dmrg_host_program_settings%nprocs > 1) dmrg_host_program_settings%runs_parallel = .true.
#else
    dmrg_host_program_settings%nprocs = 1 ! from para_info.fh
    dmrg_host_program_settings%myrank = 0 ! from para_info.fh
#endif

  end subroutine initialize_dmrg
! *********************************************************************
!  Light initialisation of DMRG parameters, used e.g. in RASSI
!  fills those (possibly) necessary fields otherwise read_dmrg_info() would use, but without reading
!  dmrg_interface.parameters
  subroutine initialize_dmrg_rassi(nstate)

    implicit none
    integer,intent(in)  :: nstate
!     integer,intent(in)  :: norb
    ! local variables
    integer              :: nproc_DMRG,lp,lcd,ierr
    character(len=20)    :: check_nr_cores
    character(len=36)    :: check_fc_offset
    character(len=1695)  :: currdir
    character(len=300)   :: project

    ! Check and set the number of CPUs
    nproc_dmrg = 1
    call getenv("QCMaquis_CPUS",check_nr_cores)
    if(check_nr_cores /= "")then
      read(check_nr_cores,'(I5)') nproc_dmrg
    else
      call getenv("OMP_NUM_THREADS",check_nr_cores)
      if(check_nr_cores /= "") read(check_nr_cores,'(I5)') nproc_dmrg
    end if

    dmrg_setup%nproc = nproc_dmrg
    !> set prefix for QCMaquis results/checkpoint files
    call getenv("Project",project)
    call getenv("CurrDir",currdir)

    dmrg_file%prefix(1:1997) = ''
    lp = 0; lcd = 0
    lp                       = len_trim(project)
    lcd                      = len_trim(currdir)
    if(trim(project) /= '' .and. trim(currdir) /= '')then
      dmrg_file%prefix(1:lp+lcd+2) = trim(currdir)//'/'//trim(project)//'.'
    else
      dmrg_file%prefix(1:2) = './'
    end if

    !> set offset for QCMaquis result/checkpoint file counting (useful if we run several RASSCF jobs in a single shot)
    dmrg_file%offset = 0
    call getenv("QCMFCO",check_fc_offset)
    if(check_fc_offset /= "")then
      read(check_fc_offset,'(I5)') dmrg_file%offset
    end if

    dmrg_state%maxroot = nstate

!   allocate(dmrg_file%qcmaquis_checkpoint_file(nstate), stat=ierr); if( ierr /= 0 ) &
!     stop ' Error in allocation: qcmaquis_checkpoint_file(:)'
!   dmrg_file%qcmaquis_checkpoint_file = ''

  end subroutine initialize_dmrg_rassi

! *********************************************************************
  subroutine print_dmrg_info(lupri,fmt2,switch,start_guess,nroots,thre)

    integer,            intent(in)    :: lupri
    integer,            intent(in)    :: switch
    integer,            intent(in)    :: nroots
    double precision,   intent(in)    :: thre
    character(len=8),   intent(in)    :: fmt2
    character(len=100), intent(inout) :: start_guess

    character(len=500)                :: mstates
    character(len=500)                :: sweeps
    character(len=500)                :: sweeps_tolerance
    character(len=500)                :: jcd_tolerance
    character(len=500)                :: jcd_maxiter
    character(len=500)                :: svd_tolerance_initial
    character(len=500)                :: svd_tolerance_final
    character(len=500)                :: orbital_ordering
    character(len=500)                :: line
    character(len=5)                  :: state_tag
    character(len=5)                  :: full_state_tag
    integer                           :: i, irootm1

    if(dmrg_host_program_settings%myrank == 0)then
      mstates               = '0'
      sweeps                = '0'
      jcd_maxiter           = '10'
      svd_tolerance_initial = '1e-50'
      orbital_ordering      = 'ascending in numerical order (default)'
      svd_tolerance_final   = ' '
      sweeps_tolerance      = ' '
      jcd_tolerance         = ' '

      write(      sweeps_tolerance,'(e9.3)') thre
      write(         jcd_tolerance,'(e9.3)') thre*0.001  ! same as molcas for Davidson
      write(   svd_tolerance_final,'(e9.3)') thre*0.001  !  in order to match Davidson

      do i = 1, size(dmrg_input%qcmaquis_input),2
        line(1:500) = dmrg_input%qcmaquis_input(i)(1:500)
        call lower_to_upper(line)
        if(trim(line) == 'TRUNCATION_INITIAL')then
          svd_tolerance_initial = trim(dmrg_input%qcmaquis_input(i+1))
        else if(trim(line) == 'TRUNCATION_FINAL')then
          svd_tolerance_final        = trim(dmrg_input%qcmaquis_input(i+1))
        else if(trim(line) == 'IETL_JCD_TOL')then
          jcd_tolerance         = trim(dmrg_input%qcmaquis_input(i+1))
        else if(trim(line) == 'IETL_JCD_MAXITER')then
           jcd_maxiter          = trim(dmrg_input%qcmaquis_input(i+1))
        else if(trim(line) == 'CONV_THRESH')then
           sweeps_tolerance     = trim(dmrg_input%qcmaquis_input(i+1))
        else if(trim(line) == 'MAX_BOND_DIMENSION')then
           mstates              = trim(dmrg_input%qcmaquis_input(i+1))
        else if(trim(line) == 'NSWEEPS')then
           sweeps               = trim(dmrg_input%qcmaquis_input(i+1))
        end if
      end do

      if(trim(mstates) == '0') mstates = 'dynamically changing (according to sweep_bond_dimensions)'

      write(lupri,fmt2//'a,t45,5x,a)') 'Number of renormalized states           ', trim(mstates)
      if(dmrg_warmup%doCIDEAS)then
        write(lupri,fmt2//'a,t45,5x,a)') 'Start guess in warm-up sweep            ','CI-DEAS'
      else
        write(lupri,fmt2//'a,t45,5x,a)') 'Start guess in warm-up sweep            ', trim(start_guess)
      end if
      write(lupri,fmt2//'a,t45,5x,a)') '(Max) number of sweeps                  ', trim(sweeps)
      write(lupri,fmt2//'a,t45,5x,a)') 'Convergence threshold (sweep tolerance) ', trim(sweeps_tolerance)
      write(lupri,fmt2//'a,t45,5x,a)') 'Jacobi-Davidson threshold               ', trim(jcd_tolerance)
      write(lupri,fmt2//'a,t45,5x,a)') 'SVD truncation threshold (initial)      ', trim(svd_tolerance_initial)
      write(lupri,fmt2//'a,t45,5x,a)') 'SVD truncation threshold (final)        ', trim(svd_tolerance_final)

      select case(switch)
        !> output before optimization
        case(1)
        if(dmrg_warmup%doFIEDLER)&
      write(lupri,fmt2//'a,t45     )') 'Fiedler vector for orbital ordering     '
        !> output after optimization
        case(2)
        do i = 1, nroots
          state_tag(1:5)      = ' '
          full_state_tag(1:5) = ' '
          irootm1        = i-1
          if(irootm1 < 10)then
            write(state_tag,'(i1)') irootm1
          else if(irootm1 < 100)then
            write(state_tag,'(i2)') irootm1
          else if(irootm1 < 1000)then
            write(state_tag,'(i3)') irootm1
          else if(irootm1 < 10000)then
            write(state_tag,'(i4)') irootm1
          else if(irootm1 < 100000)then
            write(state_tag,'(i5)') irootm1
          end if
          if(i < 10)then
            write(full_state_tag,'(i1)') i
          else if(i < 100)then
            write(full_state_tag,'(i2)') i
          else if(i < 1000)then
            write(full_state_tag,'(i3)') i
          else if(i < 10000)then
            write(full_state_tag,'(i4)') i
          else if(i < 100000)then
            write(full_state_tag,'(i5)') i
          end if
          open(898,file='internal-orbital-ordering.state.'//trim(state_tag),status='old',&
               form='formatted',action='read',position='rewind')
          orbital_ordering(1:500) = ' '
          read(898,'(a)') orbital_ordering
          write(lupri,fmt2//'a,t45,5x,a)') 'Internal orbital ordering for state   ',full_state_tag
          write(lupri,fmt2//'a,t45,5x,a  )') '                                      ',trim(orbital_ordering)
        end do
      end select
    else
      write(lupri,fmt2//'a,t45)') 'I am not master - no DMRG info print  '
    endif

  end subroutine print_dmrg_info
! *********************************************************************
  subroutine dump_dmrg_info()

    integer                         :: nproc_tmp, norb_tmp, maxroot_tmp, maxroot_save, i
    logical                         :: isthere
    real*8, allocatable             :: state_specific_energies(:)
    character(len=600), allocatable :: qcm_parameter_file(:)

    inquire(file="dmrg_interface.parameters", exist=isthere)

    if(isthere)then
      open(unit=100, status="old"    , file="dmrg_interface.parameters",         &
           action='readwrite',position='rewind',form='unformatted')
      read(100) nproc_tmp
      read(100) norb_tmp

      !> some sanity checks to make sure the same active space, orbital ordering
      !etc was used...
!should't be committed || need passby if LRras2 is activated
      if(norb_tmp /= dmrg_external%norb)                 stop 'active space  mismatch in DMRG runs'

      dmrg_setup%nproc = max(dmrg_setup%nproc,nproc_tmp)

      read(100) maxroot_tmp
      allocate(state_specific_energies(dmrg_state%maxroot+maxroot_tmp))
      state_specific_energies = 0
      read(100) state_specific_energies(1:maxroot_tmp)

      state_specific_energies(maxroot_tmp+1:dmrg_state%maxroot+maxroot_tmp) =    &
      dmrg_energy%dmrg_state_specific(1:dmrg_state%maxroot)

      !> reset the number of states and save all energies in the resized array
      maxroot_save       = dmrg_state%maxroot
      dmrg_state%maxroot = dmrg_state%maxroot + maxroot_tmp
      deallocate(dmrg_energy%dmrg_state_specific)
      allocate(dmrg_energy%dmrg_state_specific(dmrg_state%maxroot))

      dmrg_energy%dmrg_state_specific(1:dmrg_state%maxroot) =                    &
              state_specific_energies(1:dmrg_state%maxroot)

      deallocate(state_specific_energies)

      !> reset maxroot
      dmrg_state%maxroot = maxroot_save

      allocate(qcm_parameter_file(dmrg_state%maxroot+maxroot_tmp))
      qcm_parameter_file = ''
      read(100) qcm_parameter_file(1:maxroot_tmp)

      qcm_parameter_file(maxroot_tmp+1:dmrg_state%maxroot+maxroot_tmp) =    &
      dmrg_file%qcmaquis_parameter_file(1:dmrg_state%maxroot)

      !> reset the number of states and save all parameter files in the resized array
      dmrg_state%maxroot = dmrg_state%maxroot + maxroot_tmp
      deallocate(dmrg_file%qcmaquis_parameter_file)
      allocate(dmrg_file%qcmaquis_parameter_file(dmrg_state%maxroot))
      dmrg_file%qcmaquis_parameter_file(dmrg_state%maxroot) = ''

      dmrg_file%qcmaquis_parameter_file(1:dmrg_state%maxroot) =                    &
                     qcm_parameter_file(1:dmrg_state%maxroot)

      deallocate(qcm_parameter_file)

      rewind(100)

    else
      open(unit=100, status="replace", file="dmrg_interface.parameters",         &
           action='readwrite',position='rewind',form='unformatted')
    end if

    write(100) dmrg_setup%nproc
    write(100) dmrg_external%norb
    write(100) dmrg_state%maxroot
    write(100) dmrg_energy%dmrg_state_specific(1:dmrg_state%maxroot)
    write(100) dmrg_file%qcmaquis_parameter_file(1:dmrg_state%maxroot)
    write(100) dmrg_file%prefix
    close(100, status="keep")

!#ifdef _DMRG_DEBUG_
    print *, ' WRITING DMRG interface parameters'
    print *, ' ---------------------------------'
    print *, ' offset counter --> ',dmrg_file%offset
    print *, ' list of parameter files total # --> ',dmrg_state%maxroot
    do i = 1, dmrg_state%maxroot
      print *, ' parameter file for state ',i, ' --> ', trim(dmrg_file%qcmaquis_parameter_file(i))
    end do
!#endif

  end subroutine dump_dmrg_info
! *********************************************************************

  subroutine read_dmrg_info()

    integer :: ierr, i, lp, lcd
    character(len=1695)  :: UseQCMPrefix
    character(len=1695)  :: currdir
    character(len=300)   :: project

    doDMRG = .true.

    open(unit=100, status="old", file="dmrg_interface.parameters",         &
         action='readwrite',position='rewind',form='unformatted')

    read(100) dmrg_setup%nproc
    read(100) dmrg_external%norb
    read(100) dmrg_external%maxroot

    allocate(dmrg_external%dmrg_state_specific(dmrg_external%maxroot), stat=ierr); if( ierr /= 0 ) &
    stop ' Error in allocation: dmrg_state_specific(:)'
    dmrg_external%dmrg_state_specific = 0

    allocate(dmrg_file%qcmaquis_parameter_file(dmrg_external%maxroot), stat=ierr); if( ierr /= 0 ) &
    stop ' Error in allocation: qcmaquis_parameter_file(:)'
    dmrg_file%qcmaquis_parameter_file = ''

    read(100) dmrg_external%dmrg_state_specific(1:dmrg_external%maxroot)
    read(100) dmrg_file%qcmaquis_parameter_file(1:dmrg_external%maxroot)
    read(100) dmrg_file%prefix

    !> use prefix for QCMaquis results/checkpoint files only if requested
    UseQCMPrefix = ''
    call getenv("UseQCMPrefix",UseQCMPrefix)

    if(trim(UseQCMPrefix) == '')then
      dmrg_file%prefix(1:1997) = ''
      call getenv("Project",project)
      call getenv("CurrDir",currdir)
      lp = 0; lcd = 0
      lp                       = len_trim(project)
      lcd                      = len_trim(currdir)
      if(trim(project) /= '' .and. trim(currdir) /= '')then
        dmrg_file%prefix(1:lp+lcd+2) = trim(currdir)//'/'//trim(project)//'.'
      else
        dmrg_file%prefix(1:2) = './'
      end if
    end if

!#ifdef _DMRG_DEBUG_
    print *, ' READING DMRG interface parameters'
    print *, ' ---------------------------------'
    print *, ' offset counter --> ',dmrg_file%offset
    print *, ' list of parameter files total # --> ',dmrg_external%maxroot
    do i = 1, dmrg_state%maxroot
      print *, ' parameter file for state ',i, ' --> ', trim(dmrg_file%qcmaquis_parameter_file(i))
    end do
!#endif

    close(100, status="keep")

  end subroutine read_dmrg_info
! *********************************************************************

! ======================================================================
!               Save the dmrg parameters for MCLR part in MOLCAS
! ----------------------------------------------------------------------
! Input  :  maquis_model
!        :        nstate        total states that will be saved
!        :        istate        current state
!        :          nele        total electrons
!        :           MS2        Spin parameter
!        :          NRG2        Active space in DMRG
!        :        LRras2        (Reduced) Active space for CI in MCLR
!        :        result        maquis_name_results(iroot)   --  string
!        :    len_result                         --  length (integer)
! Output :  <dmrg_for_mclr.parameters>    (In scratch folder)
! ======================================================================

  subroutine save_dmrg_parameter_for_mclr(maquis_model,nstates,istate,nele,MS2,   &
                                          nrg2,lrras2,len_result,result,energy)

    character(len=3)          :: maquis_model
    integer                   :: nstates,istate,len_result,nele,MS2
    integer                   :: nrg2(20),lrras2(20)
    character(len=len_result) :: result
    real*8                    :: energy

    integer                   :: i

    if(istate.eq.1)then
      open(unit=100,file="dmrg_for_mclr.parameters")
      write(100,"(11X,L,4X,A)").true., "% doDMRG parapmters"
      write(100,"(4X,I8,4X,A)")nele,"% dmrg_state%nactel parapmters"
      write(100,"(4X,I8,4X,A)")ms2 ,"% dmrg_state%ms2    parapmters"
      do i=1,8
        write(100,"(4X,I3)",advance='no')nrg2(i)
      end do
      write(100,*)
      do i=1,8
        write(100,"(4X,I3)",advance='no')LRras2(i)
      end do
      write(100,*)
      write(100,"(4X,I8,4X,A)")nstates, "% number of calculated states"
      write(100,"(4X,A)")result(1:len_result)
      write(100,"(G20.12)")energy
    else
      write(100,"(4X,A)")result(1:len_result)
      write(100,"(G20.12)")energy
    end if

    if(istate.eq.nstates)then
      close(100)
    end if

  end subroutine save_dmrg_parameter_for_mclr

  subroutine finalize_dmrg()

    if(dmrg_host_program_settings%myrank == 0)then
      call system("rm *.pyc")
      call system("rm dmrg-input")
      ! call system("rm oneparticle.* twoparticle.* template-dmrg*")
    end if

    if(allocated(dmrg_state%iroot))                   deallocate(dmrg_state%iroot)
    if(allocated(dmrg_state%weight))                  deallocate(dmrg_state%weight)
    if(allocated(dmrg_orbital_space%initial_occ))     deallocate(dmrg_orbital_space%initial_occ)
    if(allocated(dmrg_energy%dmrg_state_specific))    deallocate(dmrg_energy%dmrg_state_specific)
    if(allocated(dmrg_energy%num_sweeps))             deallocate(dmrg_energy%num_sweeps)
    if(allocated(dmrg_energy%num_sweeps_old))         deallocate(dmrg_energy%num_sweeps_old)
    if(allocated(dmrg_energy%max_truncW))             deallocate(dmrg_energy%max_truncW)
    if(allocated(dmrg_energy%max_truncW_old))         deallocate(dmrg_energy%max_truncW_old)  
    if(allocated(dmrg_input%qcmaquis_input))          deallocate(dmrg_input%qcmaquis_input)
    if(allocated(dmrg_file%qcmaquis_parameter_file))  deallocate(dmrg_file%qcmaquis_parameter_file)
    if(allocated(dmrg_file%qcmaquis_checkpoint_file)) deallocate(dmrg_file%qcmaquis_checkpoint_file)
    if(allocated(dmrg_orbital_ordering%fiedler_order))deallocate(dmrg_orbital_ordering%fiedler_order)

  end subroutine finalize_dmrg

! ======================================================================
!       The setup of running environment for DMRG calculation
! ----------------------------------------------------------------------
! Input  : nproc
! Output : <dmrgrc> file in scratch
! ======================================================================
   subroutine set_dmrg_runtime_environment(nproc)

        integer, intent(in) :: nproc

          open(unit=121,file="dmrgrc")
            write(121,*)"--executable=""dmrg"""
            if(nproc.lt.10000.and.nproc.ge.1000)then
              write(121,"(A28,I4,A2)")"--launcher=""OMP_NUM_THREADS=",nproc ,""""
            end if
            if(nproc.lt.1000.and.nproc.ge.100)then
              write(121,"(A28,I3,A2)")"--launcher=""OMP_NUM_THREADS=",nproc ,""""
            end if
            if(nproc.lt.100.and.nproc.ge.10)then
              write(121,"(A28,I2,A2)")"--launcher=""OMP_NUM_THREADS=",nproc ,""""
            end if
            if(nproc.lt.10.and.nproc.ge.1)then
              write(121,"(A28,I1,A2)")"--launcher=""OMP_NUM_THREADS=",nproc ,""""
            end if
            write(121,*)"--get=""dmrg*"""
          close(121)

    end subroutine set_dmrg_runtime_environment

end module qcmaquis_interface_environment

! *********************************************************************
