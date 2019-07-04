!!  dmrg-interface-utils: interface to the Maquis DMRG program for various
!!                        quantum-chemistry program packages.
!!  Copyright 2013-2018 Leon Freitag, Erik Hedegaard, Sebastian Keller,
!!                      Stefan Knecht, Yingjin Ma, Christopher Stein
!!                      and Markus Reiher
!!                      Laboratory for Physical Chemistry, ETH Zurich
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

module qcmaquis_interface_cfg

! stefan: DMRG interface variables

  implicit none

  !> These are used for DMRG-SCF part
  logical, public :: doDMRG                     = .false.
  logical, public :: doVERBOSE                  = .false.
  logical, public :: doSRDFTDMRG                = .false.
  logical, public :: domcpdftDMRG               = .false.

  type qcm_warmup
    logical       :: doCIDEAS                   = .false.
    logical       :: doFiedler                  = .false.
  end type qcm_warmup
  type (qcm_warmup), save, public :: dmrg_warmup

  type qcm_orb_ordering
       character(len=1000) , allocatable :: fiedler_order(:)
  end type qcm_orb_ordering
  type (qcm_orb_ordering), save :: dmrg_orbital_ordering

  !> flag used in MPS-SI to indicate reading QCMaquis checkpoint names from rasscf.h5 files
  logical, public :: doMPSSICheckpoints         = .false.

  !> Threshold for QCMaquis, should be transfered from parent code (e.g. Molcas)
  double precision :: E_threshold               =  0.0d0

  type type_host_settings
       logical :: runs_parallel                =  .false.
       integer :: myrank                       =  0 ! rank of MPI process in host program
       integer :: nprocs                       =  1 ! number of MPI processes in host program
       character(len=7) :: dmrg_host_program   =  'molcas '
  end type type_host_settings
  type (type_host_settings), save, public :: dmrg_host_program_settings

  !> definition of "symmetry" data type
  type type_symmetry
       integer :: nirrep                        = 0
       integer :: multiplication_table(1:8,1:8) = 0
  end type type_symmetry
  type (type_symmetry), save :: dmrg_symmetry

  !> definition of DMRG parameter/input variables

  type type_setup
       character(len=2)                :: dmrg_init         = "  "
       integer                         :: nproc             =      1
  end type type_setup
  type (type_setup), save, public :: dmrg_setup

  type type_input
       character(len=500), allocatable :: qcmaquis_input(:)
       integer                         :: nr_qcmaquis_input_lines = -1
  end type type_input
  type (type_input), save, public :: dmrg_input

  type type_dmrgfiles
       character(len=1997)              :: prefix
       character(len=600), allocatable  :: qcmaquis_parameter_file(:)
       character(len=256), allocatable  :: qcmaquis_checkpoint_file(:)
       integer                          :: offset = 0
  end type type_dmrgfiles
  type (type_dmrgfiles), save, public :: dmrg_file


  !> definition of "state" data type
  type type_state
       integer              :: irefsm        = 0
       integer              :: nactel        = 0
       integer              :: mults         = 0
       integer              :: ms2           = 0
       integer              :: nroot         = 0
       integer              :: maxroot       = 0
       integer, allocatable :: iroot(:)
       real*8 , allocatable :: weight(:)
  end type type_state
  type (type_state), save :: dmrg_state

  !> definition of "orbital_space" data type
  type type_orbital_space
       integer              ::  nash(1:20)        = 0
       integer              :: LRash(1:20)        = 0        ! reduced AS for gradients (MCLR in Molcas)
       integer, allocatable :: initial_occ(:,:) ! in order to get the starting determinant for each state (Maquis)
  end type type_orbital_space
  type (type_orbital_space), save :: dmrg_orbital_space

  !> definition of "energy" data type
  type type_energy
       real*8 :: rdm                   = 0.0d0
       real*8 :: dmrg                  = 0.0d0
       real*8 , allocatable :: dmrg_state_specific(:)
       real*8 , allocatable :: max_truncW(:)
       real*8 , allocatable :: max_truncW_old(:)
       integer, allocatable :: num_sweeps(:)
       integer, allocatable :: num_sweeps_old(:)
  end type type_energy
  type (type_energy), save :: dmrg_energy

  !> DMRG-RASSI parameters
  type external_PARAMETER
       integer              :: nalpha                   = 0       ! number of alpha electrons
       integer              :: nbeta                    = 0       ! number of beta electrons
       integer              :: norb                     = 0       ! number of active orbitals
       integer              :: irrep                    = 0       ! spatial irrep
       integer              :: M                        = 0       ! number of renormalized states
       integer              :: maxroot                  = 0       ! number of states in RASSCF run
       real*8 , allocatable :: dmrg_state_specific(:)             ! energy of each state
       character(len=300)   :: masorb                   = ""      ! orbital reordering string
       logical              :: MPSrotated               = .false. ! MPSs of JOB1 and JOB2 were rotated
  end type external_PARAMETER
  type (external_PARAMETER), save, public :: dmrg_external

end module qcmaquis_interface_cfg

