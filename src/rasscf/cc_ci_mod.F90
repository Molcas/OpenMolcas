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
! Copyright (C) 2020, Oskar Weser                                      *
!***********************************************************************

#include "macros.fh"
#include "compiler_features.h"

module CC_CI_mod

use Para_Info, only: MyRank
use filesystem, only: get_errno_, getcwd_, real_path, strerror_
use linalg_mod, only: abort_, verify_
use rasscf_global, only: EMY, Ener, iter, KSDFT, lRoots, nAc, nAcPar, nAcPr2, nroots, S
use general_data, only: iSpin, nActEl, nAsh, nConf, nSym, ntot, ntot1, ntot2
use gas_data, only: iDoGas, ngssh
use generic_CI, only: CI_solver_t
use index_symmetry, only: one_el_idx, two_el_idx_flatten
use rctfld_module, only: lRF
use CI_solver_util, only: CleanMat, inv_triang_number, RDM_to_runfile, wait_and_read, write_RDM
#ifdef _ADDITIONAL_RUNTIME_CHECK_
use CI_solver_util, only: triangular_number
#endif
#ifdef _MOLCAS_MPP_
use mpi, only: MPI_COMM_WORLD, MPI_REAL8
use Para_Info, only: Is_Real_Par
use Definitions, only: MPIInt
#endif
use Constants, only: Zero, Two, Half
use Definitions, only: wp, iwp, u6

implicit none
private

public :: Do_CC_CI, CC_CI_solver_t, write_RDM

logical(kind=iwp) :: Do_CC_CI = .false.

interface
  function isfreeunit(iseed)
    import :: iwp
    integer(kind=iwp) :: isfreeunit
    integer(kind=iwp), intent(in) :: iseed
  end function
end interface

type, extends(CI_solver_t) :: CC_CI_solver_t
contains
  procedure :: run => CC_CI_ctl
  procedure :: cleanup
end type

interface CC_CI_solver_t
  module procedure :: construct_CC_CI_solver_t
end interface

contains

subroutine CC_CI_ctl(this,actual_iter,ifinal,iroot,weight,CMO,DIAF,D1I_AO,D1A_AO,TUVX,F_IN,D1S_MO,DMAT,PSMAT,PAMAT)

  use fcidump_reorder, only: get_P_GAS, get_P_inp, ReOrFlag, ReOrInp
  use fcidump, only: make_fcidumps, transform

  class(CC_CI_solver_t), intent(in) :: this
  integer, intent(in) :: actual_iter, ifinal, iroot(nroots)
  real(kind=wp), intent(in) :: weight(nroots), CMO(nTot2), DIAF(nTot), D1I_AO(nTot2), D1A_AO(nTot2), TUVX(nAcpr2)
  real(kind=wp), intent(inout) :: F_In(nTot1), D1S_MO(nAcPar)
  real(kind=wp), intent(out) :: DMAT(nAcpar), PSMAT(nAcpr2), PAMAT(nAcpr2)
  integer(kind=iwp) :: jRoot
# ifdef _MOLCAS_MPP_
  integer(kind=MPIInt) :: error
# endif
  real(kind=wp) :: energy(nroots), folded_Fock(nAcPar), orbital_E(nTot)
  integer(kind=iwp), allocatable :: permutation(:)
  character(len=*), parameter :: ascii_fcidmp = 'FCIDUMP', h5_fcidmp = 'H5FCIDUMP'

  unused_var(this)
  unused_var(ifinal)
  unused_var(weight)

  if (size(iroot) >= 2) call abort_('SA-CC-CASSCF yet to be implemented.')

  ! SOME DIRTY SETUPS
  S = Half*real(iSpin-1,kind=wp)

  call check_options(lRoots,lRf,KSDFT,iDoGAS)

  ! Produce a working FCIDUMP file
  if (ReOrFlag /= 0) then
    allocate(permutation(sum(nAsh(:nSym))))
    if (ReOrFlag >= 2) permutation(:) = get_P_inp(ReOrInp)
    if (ReOrFlag == -1) permutation(:) = get_P_GAS(nGSSH)
  end if

  ! This call is not side effect free, sets EMY and modifies F_IN
  call transform(actual_iter,CMO,DIAF,D1I_AO,D1A_AO,D1S_MO,F_IN,orbital_E,folded_Fock)

  ! Fortran Standard 2008 12.5.2.12:
  ! Allocatable actual arguments that are passed to
  ! non-allocatable, optional dummy arguments are **not** present.
  call make_fcidumps(ascii_fcidmp,h5_fcidmp,orbital_E,folded_Fock,TUVX,EMY,permutation)

  ! Run CC
# ifdef _MOLCAS_MPP_
  if (is_real_par()) call MPI_Barrier(MPI_COMM_WORLD,error)
# endif

  call run_CC_CI(ascii_fcidmp,h5_fcidmp,fake_run=actual_iter == 1,energy=energy,D1S_MO=D1S_MO,DMAT=DMAT,PSMAT=PSMAT,PAMAT=PAMAT)
  do jRoot=1,lRoots
    ENER(jRoot,ITER) = energy(jRoot)
  end do

  if (nAsh(1) /= nac) call dblock(dmat)

end subroutine CC_CI_ctl

subroutine run_CC_CI(ascii_fcidmp,h5_fcidmp,fake_run,energy,D1S_MO,DMAT,PSMAT,PAMAT)

  character(len=*), intent(in) :: ascii_fcidmp, h5_fcidmp
  logical(kind=iwp), intent(in) :: fake_run
  real(kind=wp), intent(out) :: energy(nroots), D1S_MO(nAcPar), DMAT(nAcpar), PSMAT(nAcpr2), PAMAT(nAcpr2)
  real(kind=wp) :: previous_energy(nroots)
  character(len=*), parameter :: energy_file = 'NEWCYCLE', input_name = 'CC_CI.inp'

  if (fake_run) then
    energy(:) = previous_energy(:)
  else
    call make_inp(input_name)
    if (myrank == 0) call write_user_message(input_name,ascii_fcidmp,h5_fcidmp)
    call wait_and_read(energy_file,energy)
    previous_energy = energy
  end if
  call read_CC_RDM(DMAT,D1S_MO,PSMAT,PAMAT)
  call RDM_to_runfile(DMAT,D1S_MO,PSMAT,PAMAT)

end subroutine run_CC_CI

subroutine make_inp(input_name)

  character(len=*), intent(in) :: input_name

  write(u6,*) input_name
  call abort_('make_inp has to be implemented.')

end subroutine

subroutine cleanup(this)

  use fcidump, only: fcidump_cleanup => cleanup

  class(CC_CI_solver_t), intent(inout) :: this

  unused_var(this)
  call fcidump_cleanup()

end subroutine cleanup

function construct_CC_CI_solver_t() result(res)

  type(CC_CI_solver_t) :: res

  unused_var(res)
# if defined (_WARNING_WORKAROUND_) && ! defined (EMPTY_TYPE_INIT)
  res = res
# endif
  ! Due to possible size of active space, arrays of nConf size
  ! need to be avoided.  For this reason set nConf to zero.
  write(u6,*) ' DCC-CI activated. List of Confs might get lengthy.'
  write(u6,*) ' Number of Configurations computed by GUGA: ',nConf
  write(u6,*) ' nConf variable is set to zero to avoid JOBIPH i/o'
  nConf = 0

end function

subroutine check_options(lroots,lRf,KSDFT,DoGAS)

  integer(kind=iwp), intent(in) :: lroots
  logical(kind=iwp), intent(in) :: lRf, DoGAS
  character(len=*), intent(in) :: KSDFT
  logical(kind=iwp) :: Do_ESPF

  call verify_(lroots == 1,"CC-CI doesn't support State Average!")

  call DecideOnESPF(Do_ESPF)
  if (lRf .or. (KSDFT /= 'SCF') .or. Do_ESPF) call abort_('CC CI does not support Reaction Field yet!')

  call verify_(.not. DoGAS,'CC CI does not support GASSCF yet!')

end subroutine check_options

subroutine write_user_message(input_name,ascii_fcidmp,h5_fcidmp)

  character(len=*), intent(in) :: input_name, ascii_fcidmp, h5_fcidmp
  integer(kind=iwp) :: err
  character(len=1024) :: WorkDir

  call getcwd_(WorkDir,err)
  if (err /= 0) write(u6,*) strerror_(get_errno_())

  write(u6,'(A)') 'Run coupled cluster CI externally.'
  write(u6,'(A)') 'Get the (example) coupled cluster input:'
  write(u6,'(4x, A, 1x, A, 1x, A)') 'cp',real_path(input_name),'$CC_RUN_DIR'
  write(u6,'(A)') 'Get the ASCII formatted FCIDUMP:'
  write(u6,'(4x, A, 1x, A, 1x, A)') 'cp',real_path(ascii_fcidmp),'$CC_RUN_DIR'
  write(u6,'(A)') 'Or the HDF5 FCIDUMP:'
  write(u6,'(4x, A, 1x, A, 1x, A)') 'cp',real_path(h5_fcidmp),'$CC_RUN_DIR'
  write(u6,*)
  write(u6,'(A)') 'When finished do:'
  ! TODO(Oskar, Thomas): Change accordingly
  write(u6,'(4x, A)') 'cp PSMAT.dat PAMAT.dat '//trim(WorkDir)
  write(u6,'(4x, A)') 'echo $your_RDM_Energy > '//real_path('NEWCYCLE')
  call xflush(u6)

end subroutine write_user_message

!>  @brief
!>    Read DCC RDM files
!>
!>  @author Oskar Weser
!>
!>  @param[out] DMAT Average 1 body density matrix
!>  @param[out] D1S_MO Average spin 1-dens matrix
!>  @param[out] PSMAT Average symm. 2-dens matrix
!>  @param[out] PAMAT Average antisymm. 2-dens matrix
subroutine read_CC_RDM(DMAT,D1S_MO,PSMAT,PAMAT)

  real(kind=wp), intent(out) :: DMAT(nAcpar), D1S_MO(nAcPar), PSMAT(nAcpr2), PAMAT(nAcpr2)
# ifdef _MOLCAS_MPP_
  integer(kind=MPIInt) :: error
# endif

  if (myrank == 0) then
    call read_2RDM('PSMAT.dat',PSMAT)
    call read_2RDM('PAMAT.dat',PAMAT)
    call calc_1RDM(PSMAT,DMAT)
    call cleanMat(DMAT)
  end if
  ! No spin resolved RDMs available at the moment.
  D1S_MO(:) = Zero
  ! Could be changed into non blocking BCast
# ifdef _MOLCAS_MPP_
  if (is_real_par()) then
    call MPI_Bcast(PSMAT,size(PSMAT,kind=MPIInt),MPI_REAL8,0_MPIInt,MPI_COMM_WORLD,error)
    call MPI_Bcast(PAMAT,size(PAMAT,kind=MPIInt),MPI_REAL8,0_MPIInt,MPI_COMM_WORLD,error)
    call MPI_Bcast(DMAT,size(DMAT,kind=MPIInt),MPI_REAL8,0_MPIInt,MPI_COMM_WORLD,error)
  end if
# endif

end subroutine read_CC_RDM

subroutine calc_1RDM(PSMAT,DMAT)

  real(kind=wp), intent(in) :: PSMAT(:)
  real(kind=wp), intent(out) :: DMAT(:)
  integer(kind=iwp) :: p, pq, q, r
  debug_function_name('calc_1RDM')

  ASSERT(size(PSMAT) == triangular_number(size(DMAT)))

  DMAT(:) = Zero
  do pq=1,size(DMAT)
    call one_el_idx(pq,p,q)
    do r=1,inv_triang_number(size(DMAT))
      DMAT(pq) = DMAT(pq)+PSMAT(two_el_idx_flatten(p,q,r,r))
    end do
  end do
  DMAT(:) = DMAT(:)*Two/real(nActEl-1,kind=wp)

end subroutine calc_1RDM

subroutine read_2RDM(path,RDM_2)

  character(len=*), intent(in) :: path
  real(kind=wp), intent(out) :: RDM_2(:)
  integer(kind=iwp) :: curr_line, file_id, i, io_err, n_lines
  integer(kind=iwp), parameter :: arbitrary_magic_number = 42

  ASSERT(size(RDM_2) == nAcpr2)
  n_lines = inv_triang_number(nAcpr2)

  file_id = isfreeunit(arbitrary_magic_number)
  i = 1
  call molcas_open(file_id,trim(path))
  do curr_line=1,n_lines
    read(file_id,*,iostat=io_err) RDM_2(i:i+curr_line-1)
    call verify_(io_err == 0,'Error on reading 2-RDMs.')
    i = i+curr_line
  end do
  close(file_id)

end subroutine read_2RDM

end module CC_CI_mod
