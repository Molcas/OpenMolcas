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
* Copyright (C) 2020, Oskar Weser                                      *
************************************************************************
      module CC_CI_mod
#ifdef _MOLCAS_MPP_
      use mpi
#endif
#ifdef NAGFOR
      use f90_unix_proc, only : sleep
#endif
      use filesystem, only : chdir_, getcwd_, get_errno_, strerror_,
     &    real_path
      use fortran_strings, only : str
      use stdalloc, only : mma_allocate, mma_deallocate, mxMem

      use rasscf_data, only : lRoots, nRoots, iRoot
      use general_data, only : nSym, nConf

      use CI_solver_util, only: wait_and_read, abort_

      implicit none
      save
      private
      public :: CC_CI_ctl, Do_CC_CI, init, cleanup
      logical :: Do_CC_CI = .false.
#include "para_info.fh"
#ifdef _MOLCAS_MPP_
      integer*4 :: error
#endif


      interface
        integer function isfreeunit(iseed)
          integer, intent(in) :: iseed
        end function
      end interface
      contains

      subroutine CC_CI_ctl(actual_iter, CMO, DIAF, D1I_AO, D1A_AO,
     &                 TUVX, F_IN, D1S_MO, DMAT, PSMAT, PAMAT)
      use general_data, only : iSpin, ntot, ntot1, ntot2, nAsh, nBas
      use rasscf_data, only : iter, lRoots, nRoots, S, KSDFT, EMY,
     &    rotmax, Ener, Nac, nAcPar, nAcpr2

      use gugx_data, only : IfCAS
      use gas_data, only : ngssh, iDoGas, nGAS, iGSOCCX

      use fcidump_reorder, only : get_P_GAS, get_P_inp,ReOrFlag,ReOrInp
      use fcidump, only : make_fcidumps, transform

#include "output_ras.fh"
#include "rctfld.fh"
#include "timers.fh"
      integer, intent(in) :: actual_iter
      real*8, intent(in) ::
     &    CMO(nTot2), DIAF(nTot),
     &    D1I_AO(nTot2), D1A_AO(nTot2), TUVX(nAcpr2)
      real*8, intent(inout) :: F_In(nTot1), D1S_MO(nAcPar)
      real*8, intent(out) :: DMAT(nAcpar),
     &    PSMAT(nAcpr2), PAMAT(nAcpr2)
      real*8, save :: energy
      integer :: iPRLEV, iOff, iSym, iBas, i, j, jRoot
      integer, allocatable :: permutation(:)
      real*8 :: orbital_E(nTot), folded_Fock(nAcPar)

      parameter(ROUTINE = 'CC_CI_ctl')
      character(*), parameter ::
     &  ascii_fcidmp = 'FCIDUMP', h5_fcidmp = 'H5FCIDUMP'

      call qEnter(routine)

! SOME DIRTY SETUPS
      S = 0.5d0 * dble(iSpin - 1)

      call check_options(lRoots, lRf, KSDFT, iDoGAS)

! Produce a working FCIDUMP file
      if (ReOrFlag /= 0) then
        allocate(permutation(sum(nAsh(:nSym))))
        if (ReOrFlag >= 2) permutation(:) = get_P_inp(ReOrInp)
        if (ReOrFlag == -1) permutation(:) = get_P_GAS(nGSSH)
      end if

! This call is not side effect free, sets EMY and modifies F_IN
      call transform(iter, CMO, DIAF, D1I_AO, D1A_AO, D1S_MO,
     &      F_IN, orbital_E, folded_Fock)

! Fortran Standard 2008 12.5.2.12:
! Allocatable actual arguments that are passed to
! non-allocatable, optional dummy arguments are **not** present.
      call make_fcidumps(ascii_fcidmp, h5_fcidmp,
     &                   orbital_E, folded_Fock, TUVX, EMY, permutation)

! Run NECI
      call Timing(Rado_1, Swatch, Swatch, Swatch)
#ifdef _MOLCAS_MPP_
      if (is_real_par()) call MPI_Barrier(MPI_COMM_WORLD, error)
#endif

      call run_CC_CI(ascii_fcidmp, h5_fcidmp,
     &      fake_run=actual_iter == 1, energy=energy,
     &      D1S_MO=D1S_MO, DMAT=DMAT, PSMAT=PSMAT, PAMAT=PAMAT)
! NECIen so far is only the energy for the GS.
! Next step it will be an array containing energies for all the optimized states.
      do jRoot = 1, lRoots
        ENER(jRoot, ITER) = energy
      end do

      if (nAsh(1) /= nac) call dblock(dmat)


      call Timing(Rado_2, Swatch, Swatch, Swatch)
      Rado_2 = Rado_2 - Rado_1
      Rado_3 = Rado_3 + Rado_2

      call qExit(routine)
      end subroutine CC_CI_ctl


      subroutine run_CC_CI(ascii_fcidmp, h5_fcidmp,
     &      fake_run, energy, D1S_MO, DMAT, PSMAT, PAMAT)
        use rasscf_data, only : nAcPar, nAcPr2
        implicit none
        character(*), intent(in) :: ascii_fcidmp, h5_fcidmp
        logical, intent(in) :: fake_run
        real*8, intent(out) :: energy, D1S_MO(nAcPar), DMAT(nAcpar),
     &      PSMAT(nAcpr2), PAMAT(nAcpr2)
        real*8, save :: previous_energy = 0.0d0

        character(*), parameter :: input_name = 'CC_CI.inp'

        if (fake_run) then
          energy = previous_energy
        else
          call make_inp(input_name)
          if (myrank == 0) then
            call write_user_message(input_name, ascii_fcidmp, h5_fcidmp)
          end if
          call wait_and_read(energy)
          previous_energy = energy
        end if
        call get_CC_RDM(D1S_MO, DMAT, PSMAT, PAMAT)
      end subroutine run_CC_CI

      subroutine make_inp(input_name)
        character(*), intent(in) :: input_name
        call abort_('make_inp has to be implemented.')
      end subroutine


      subroutine cleanup()
        use fcidump, only : fcidump_cleanup => cleanup
        call fcidump_cleanup()
      end subroutine cleanup

      subroutine init()
! Due to possible size of active space arrays of nConf
! size need to be avoided.  For this reason set nConf to zero.
        write(6,*) ' NECI activated. List of Confs might get lengthy.'
        write(6,*) ' Number of Configurations computed by GUGA: ', nConf
        write(6,*) ' nConf variable is set to zero to avoid JOBIPH i/o'
        nConf= 0
      end subroutine


      subroutine check_options(lroots, lRf, KSDFT, DoGAS)
        integer, intent(in) :: lroots
        logical, intent(in) :: lRf, DoGAS
        character(*), intent(in) :: KSDFT
        logical :: Do_ESPF
        if (lroots > 1) then
          call abort_('CC CI does not support State Average yet!')
        end if
        call DecideOnESPF(Do_ESPF)
        if ( lRf .or. KSDFT /= 'SCF' .or. Do_ESPF) then
          call abort_('CC CI does not support Reaction Field yet!')
        end if
        if (DoGAS) call abort_('GAS not yet supported!')
      end subroutine check_options

      subroutine write_user_message(
     &      input_name, ascii_fcidmp, h5_fcidmp)
        character(*), intent(in) :: input_name, ascii_fcidmp, h5_fcidmp
        character(1024) :: WorkDir
        integer :: err

        call getcwd_(WorkDir, err)
        if (err /= 0) write(6, *) strerror_(get_errno_())

        write(6,'(A)')'Run coupled cluster CI externally.'
        write(6,'(A)')'Get the (example) coupled cluster input:'
        write(6,'(4x, A, 1x, A, 1x, A)')
     &    'cp', real_path(input_name), '$NECI_RUN_DIR'
        write(6,'(A)')'Get the ASCII formatted FCIDUMP:'
        write(6,'(4x, A, 1x, A, 1x, A)')
     &    'cp', real_path(ascii_fcidmp), '$NECI_RUN_DIR'
        write(6,'(A)')'Or the HDF5 FCIDUMP:'
        write(6,'(4x, A, 1x, A, 1x, A)')
     &    'cp', real_path(h5_fcidmp), '$NECI_RUN_DIR'
        write(6, *)
        write(6,'(A)') "When finished do:"
! TODO(Oskar, Thomas): Change accordingly
        write(6,'(4x, A)')
     &    'cp TwoRDM_aaaa.1 TwoRDM_abab.1 TwoRDM_abba.1 '//
     &    'TwoRDM_bbbb.1 TwoRDM_baba.1 TwoRDM_baab.1 '//trim(WorkDir)
        write(6,'(4x, A)')
     &    'echo $your_RDM_Energy > '//real_path('NEWCYCLE')
        call xflush(6)
      end subroutine write_user_message

!> Generate density matrices for Molcas
!>   Neci density matrices are stored in Files TwoRDM_**** (in spacial orbital basis).
!>   I will be reading them from those formatted files for the time being.
!>   Next it will be nice if NECI prints them out already in Molcas format.
      subroutine get_CC_RDM(D1S_MO, DMAT, PSMAT, PAMAT)
        use general_data, only : JobIPH
        use rasscf_data, only : iAdr15, Weight, nAcPar, nAcPr2
        use fciqmc_read_RDM, only : read_neci_RDM
        implicit none
        real*8, intent(out) ::
     &      D1S_MO(nAcPar), DMAT(nAcpar),
     &      PSMAT(nAcpr2), PAMAT(nAcpr2)
        real*8, allocatable ::
!> one-body density
     &    DTMP(:),
!> symmetric two-body density
     &    Ptmp(:),
!> antisymmetric two-body density
     &    PAtmp(:),
!> one-body spin density
     &    DStmp(:)
        real*8 :: Scal
        integer :: jRoot, kRoot, iDisk, jDisk

        call mma_allocate(DTMP, nAcPar, label='Dtmp ')
        call mma_allocate(DStmp, nAcPar, label='DStmp')
        call mma_allocate(Ptmp, nAcPr2, label='Ptmp ')
        call mma_allocate(PAtmp, nAcPr2, label='PAtmp')

        call read_neci_RDM(DTMP, DStmp, Ptmp, PAtmp)

! Compute average density matrices
        do jRoot = 1, lRoots
          Scal = 0.0d0
          do kRoot = 1, nRoots
            if (iRoot(kRoot) == jRoot) Scal = Weight(kRoot)
          end do
          DMAT(:) = SCAL * DTMP(:)
          D1S_MO(:) = SCAL * PSMAT(:)
          PSMAT(:) = SCAL * Ptmp(:)
          PAMAT(:) = SCAL * PAtmp(:)
! Put it on the RUNFILE
          call Put_D1MO(DTMP,NACPAR)
          call Put_P2MO(Ptmp,NACPR2)
! Save density matrices on disk
          iDisk = IADR15(4)
          jDisk = IADR15(3)
          call DDafile(JOBIPH, 1, DTMP, NACPAR, jDisk)
          call DDafile(JOBIPH, 1, DStmp, NACPAR, jDisk)
          call DDafile(JOBIPH, 1, Ptmp, NACPR2, jDisk)
          call DDafile(JOBIPH, 1, PAtmp, NACPR2, jDisk)
        end do

        call mma_deallocate(DTMP)
        call mma_deallocate(DStmp)
        call mma_deallocate(Ptmp)
        call mma_deallocate(PAtmp)
      end subroutine get_CC_RDM

      end module CC_CI_mod
