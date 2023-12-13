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
! Copyright (C) 2017, Leon Freitag                                     *
!               2017, Stefan Knecht                                    *
!***********************************************************************
!***************
! This is the wrapper for Celestino Angeli's NEVPT2 code, modified for
! DMRG and Cholesky support
! Author: Leon Freitag, ETH Zurich
! last modifications: S. Knecht, ETH Zurich, August 2017
!***************

subroutine nevpt2(iReturn)

use koopro4QD, only: koopro4qd_driver
use qdnevpt_core, only: qdnevpt
use nevpt2wfn, only: nevpt2wfn_estore                                  ! wfn file
use nevpt_header, only: print_nevpt_header                             ! header
use nevpt2_cfg, only: do_dmrg_pt, no_pc, nr_states, skip_koopro_molcas ! input variables
use info_state_energy, only: e2en, e2mp                                ! energies
use Para_Info, only: King
use Definitions, only: iwp, u6, r4

implicit none
integer(kind=iwp), intent(out) :: iReturn
logical(kind=iwp) :: nevpth5_exists = .true.
character(len=256) :: refwfnfile
! dummy variables for the koopro4QD call
real(kind=r4) :: t13
logical(kind=iwp), parameter :: rel_ham = .false., from_molcas = .true.

! Make sure we run single-threaded in an MPI environment

if (KING()) then
  !> initialize
  refwfnfile = ''

  !> print NEVPT2 header
  do_dmrg_pt = .true.
  call print_nevpt_header(rel_ham)
  call xflush(u6)

  !> read and process input
  call rdinput(refwfnfile)
  call xflush(u6)

  !> read and process input
  call pt2init(refwfnfile)
  call xflush(u6)

  if (.not. skip_koopro_molcas) then
    write(u6,'(/a )') ' Calculating the Koopmans'' matrices'
    write(u6,'( a )') ' ----------------------------------'
    call koopro4QD_driver(rel_ham,t13,from_molcas)
    call xflush(u6)
  else
    write(u6,*) 'Skipping the calculation of Koopmans'' matrices.'
    call f_inquire('nevpt.h5', nevpth5_exists)
    if (.not. nevpth5_exists) then
      call WarningMessage(1,'Requested to skip the Koopmans'' matrix calculation step but nevpt.h5 file not found!')
      call Abend()
    end if
  end if
  write(u6,'(/a )') ' Starting NEVPT2 perturbation summation'
  write(u6,'( a )') ' --------------------------------------'
  call qdnevpt(rel_ham,t13,from_molcas)
  call xflush(u6)

  !> store the PT2 energy and effective hamiltonian on the wavefunction file
  call nevpt2wfn_estore()

  ! Test verification, a tolerance of 10^-6 is assumed for the test
  call Add_Info('E_SC_NEVPT2',e2mp,nr_states,6)
  if (.not. no_pc) then
    call Add_Info('E_PC_NEVPT2',e2en,nr_states,6)
  end if

  !> clean up and free memory
  call pt2close()
end if
iReturn = 0

end subroutine nevpt2
