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

module CC_CI_mod
#ifdef _MOLCAS_MPP_
    use mpi
#endif
    use definitions, only: MPIInt, wp
    use filesystem, only: chdir_, getcwd_, get_errno_, strerror_, &
     real_path
    use linalg_mod, only: verify_, abort_
    use fortran_strings, only: str
    use stdalloc, only : mma_allocate, mma_deallocate, mxMem

    use rasscf_data, only: iter, lRoots, nRoots, iRoot, EMY, &
         S, KSDFT, rotmax, Ener, iAdr15, Weight, nAc, nAcPar, nAcPr2
    use general_data, only: iSpin, nSym, nConf, JobIPH, &
         ntot, ntot1, ntot2, nAsh, nBas, nActEl
    use gugx_data, only: IfCAS
    use gas_data, only: ngssh, iDoGas, nGAS, iGSOCCX

    use generic_CI, only: CI_solver_t
    use index_symmetry, only: one_el_idx, two_el_idx, &
        one_el_idx_flatten, two_el_idx_flatten
    use CI_solver_util, only: wait_and_read, RDM_to_runfile, &
        CleanMat

    implicit none
    save
    private
    public :: Do_CC_CI, CC_CI_solver_t, write_RDM
    logical :: Do_CC_CI = .false.

#include "para_info.fh"
#include "rctfld.fh"

    interface
        integer function isfreeunit(iseed)
            integer, intent(in) :: iseed
        end function
    end interface

    type, extends(CI_solver_t) :: CC_CI_solver_t
        contains
        procedure, nopass :: init
        procedure, nopass :: run => CC_CI_ctl
        procedure, nopass :: cleanup
    end type

contains

    subroutine CC_CI_ctl(actual_iter, CMO, DIAF, D1I_AO, D1A_AO, &
                         TUVX, F_IN, D1S_MO, DMAT, PSMAT, PAMAT)
        use fcidump_reorder, only : get_P_GAS, get_P_inp,ReOrFlag,ReOrInp
        use fcidump, only : make_fcidumps, transform
        integer, intent(in) :: actual_iter
        real(wp), intent(in) :: &
            CMO(nTot2), DIAF(nTot), D1I_AO(nTot2), D1A_AO(nTot2), TUVX(nAcpr2)
        real(wp), intent(inout) :: F_In(nTot1), D1S_MO(nAcPar)
        real(wp), intent(out) :: DMAT(nAcpar), PSMAT(nAcpr2), PAMAT(nAcpr2)
        real(wp) :: energy
        integer :: jRoot
        integer, allocatable :: permutation(:)
        real(wp) :: orbital_E(nTot), folded_Fock(nAcPar)
#ifdef _MOLCAS_MPP_
        integer(MPIInt) :: error
#endif
        character(len=*), parameter :: &
            ascii_fcidmp = 'FCIDUMP', h5_fcidmp = 'H5FCIDUMP'


! SOME DIRTY SETUPS
        S = 0.5_wp * dble(iSpin - 1)

        call check_options(lRoots, lRf, KSDFT, iDoGAS)

! Produce a working FCIDUMP file
        if (ReOrFlag /= 0) then
            allocate(permutation(sum(nAsh(:nSym))))
            if (ReOrFlag >= 2) permutation(:) = get_P_inp(ReOrInp)
            if (ReOrFlag == -1) permutation(:) = get_P_GAS(nGSSH)
        end if

! This call is not side effect free, sets EMY and modifies F_IN
        call transform(actual_iter, CMO, DIAF, D1I_AO, D1A_AO, D1S_MO, F_IN, orbital_E, folded_Fock)

! Fortran Standard 2008 12.5.2.12:
! Allocatable actual arguments that are passed to
! non-allocatable, optional dummy arguments are **not** present.
        call make_fcidumps(ascii_fcidmp, h5_fcidmp, &
                         orbital_E, folded_Fock, TUVX, EMY, permutation)

! Run CC
#ifdef _MOLCAS_MPP_
        if (is_real_par()) call MPI_Barrier(MPI_COMM_WORLD, error)
#endif

        call run_CC_CI(ascii_fcidmp, h5_fcidmp, &
                       fake_run=actual_iter == 1, energy=energy, &
                       D1S_MO=D1S_MO, DMAT=DMAT, PSMAT=PSMAT, PAMAT=PAMAT)
        do jRoot = 1, lRoots
            ENER(jRoot, ITER) = energy
        end do

        if (nAsh(1) /= nac) call dblock(dmat)

    end subroutine CC_CI_ctl


    subroutine run_CC_CI(ascii_fcidmp, h5_fcidmp, &
                         fake_run, energy, D1S_MO, DMAT, PSMAT, PAMAT)
        character(len=*), intent(in) :: ascii_fcidmp, h5_fcidmp
        logical, intent(in) :: fake_run
        real(wp), intent(out) :: energy, D1S_MO(nAcPar), DMAT(nAcpar), &
            PSMAT(nAcpr2), PAMAT(nAcpr2)
        real(wp), save :: previous_energy = 0.0_wp

        character(len=*), parameter :: &
            input_name = 'CC_CI.inp', energy_file = 'NEWCYCLE'

        if (fake_run) then
            energy = previous_energy
        else
            call make_inp(input_name)
            if (myrank == 0) call write_user_message(input_name, ascii_fcidmp, h5_fcidmp)
            call wait_and_read(energy_file, energy)
            previous_energy = energy
        end if
        call read_CC_RDM(DMAT, D1S_MO, PSMAT, PAMAT)
        call RDM_to_runfile(DMAT, D1S_MO, PSMAT, PAMAT)
    end subroutine run_CC_CI

    subroutine make_inp(input_name)
        character(len=*), intent(in) :: input_name
        write(6, *) input_name
        call abort_('make_inp has to be implemented.')
    end subroutine


    subroutine cleanup()
        use fcidump, only : fcidump_cleanup => cleanup
        call fcidump_cleanup()
    end subroutine cleanup

    subroutine init()
! Due to possible size of active space arrays of nConf
! size need to be avoided.  For this reason set nConf to zero.
        write(6,*) ' DCC-CI activated. List of Confs might get lengthy.'
        write(6,*) ' Number of Configurations computed by GUGA: ', nConf
        write(6,*) ' nConf variable is set to zero to avoid JOBIPH i/o'
        nConf= 0
    end subroutine


    subroutine check_options(lroots, lRf, KSDFT, DoGAS)
        integer, intent(in) :: lroots
        logical, intent(in) :: lRf, DoGAS
        character(len=*), intent(in) :: KSDFT
        logical :: Do_ESPF
        call verify_(lroots == 1, "CC-CI doesn't support State Average!")

        call DecideOnESPF(Do_ESPF)
        if ( lRf .or. KSDFT /= 'SCF' .or. Do_ESPF) then
            call abort_('CC CI does not support Reaction Field yet!')
        end if

        call verify_(.not. DoGAS, 'CC CI does not support GASSCF yet!')
    end subroutine check_options

    subroutine write_user_message(input_name, ascii_fcidmp, h5_fcidmp)
        character(len=*), intent(in) :: input_name, ascii_fcidmp, h5_fcidmp
        character(len=1024) :: WorkDir
        integer :: err

        call getcwd_(WorkDir, err)
        if (err /= 0) write(6, *) strerror_(get_errno_())

        write(6,'(A)')'Run coupled cluster CI externally.'
        write(6,'(A)')'Get the (example) coupled cluster input:'
        write(6,'(4x, A, 1x, A, 1x, A)') &
                'cp', real_path(input_name), '$CC_RUN_DIR'
        write(6,'(A)')'Get the ASCII formatted FCIDUMP:'
        write(6,'(4x, A, 1x, A, 1x, A)') &
                'cp', real_path(ascii_fcidmp), '$CC_RUN_DIR'
        write(6,'(A)')'Or the HDF5 FCIDUMP:'
        write(6,'(4x, A, 1x, A, 1x, A)') &
                'cp', real_path(h5_fcidmp), '$CC_RUN_DIR'
        write(6, *)
        write(6,'(A)') "When finished do:"
! TODO(Oskar, Thomas): Change accordingly
        write(6,'(4x, A)') 'cp PSMAT.dat PAMAT.dat '//trim(WorkDir)
        write(6,'(4x, A)') 'echo $your_RDM_Energy > '//real_path('NEWCYCLE')
        call xflush(6)
    end subroutine write_user_message

!>  @brief
!>    Read DCC RDM files
!>
!>  @author Oskar Weser
!>
!>  @paramin[out] DMAT Average 1 body density matrix
!>  @paramin[out] DSPN Average spin 1-dens matrix
!>  @paramin[out] PSMAT Average symm. 2-dens matrix
!>  @paramin[out] PAMAT Average antisymm. 2-dens matrix
    subroutine read_CC_RDM(DMAT, D1S_MO, PSMAT, PAMAT)
        real(wp), intent(out) :: DMAT(nAcpar), D1S_MO(nAcPar), &
                                 PSMAT(nAcpr2), PAMAT(nAcpr2)
#ifdef _MOLCAS_MPP_
        integer(MPIInt) :: error
#endif
        if (myrank == 0) then
            call read_2RDM('PSMAT.dat', PSMAT)
            call read_2RDM('PAMAT.dat', PAMAT)
            call calc_1RDM(PSMAT, DMAT)
            call cleanMat(DMAT)
        end if
! No spin resolved RDMs available at the moment.
        D1S_MO(:) = 0.0
! Could be changed into non blocking BCast
#ifdef _MOLCAS_MPP_
        if (is_real_par()) then
            call MPI_Bcast(PSMAT, size(PSMAT, kind=MPIInt), &
                           MPI_REAL8, 0_MPIInt, MPI_COMM_WORLD, error)
            call MPI_Bcast(PAMAT, size(PAMAT, kind=MPIInt), &
                           MPI_REAL8, 0_MPIInt, MPI_COMM_WORLD, error)
            call MPI_Bcast(DMAT, size(DMAT, kind=MPIInt), &
                           MPI_REAL8, 0_MPIInt, MPI_COMM_WORLD, error)
        end if
#endif
    end subroutine read_CC_RDM

    subroutine calc_1RDM(PSMAT, DMAT)
        real(wp), intent(in) :: PSMAT(:)
        real(wp), intent(out) :: DMAT(:)
        debug_function_name("calc_1RDM")

        integer :: pq, p, q, r

        ASSERT(size(PSMAT) == triangular_number(size(DMAT)))

        DMAT(:) = 0._wp
        do pq = lbound(DMAT, 1), ubound(DMAT, 1)
            call one_el_idx(pq, p, q)
            do r = 1, inv_triang_number(size(DMAT))
                DMAT(pq) = DMAT(pq) + PSMAT(two_el_idx_flatten(p, q, r, r))
            end do
        end do
        DMAT(:) = DMAT(:) * 2.0_wp / real(nActEl - 1, kind=wp)
    end subroutine

    subroutine read_2RDM(path, RDM_2)
        character(len=*), intent(in) :: path
        real(wp), intent(out) :: RDM_2(:)

        integer :: file_id, io_err, curr_line, i, n_lines
        integer, parameter :: arbitrary_magic_number = 42


        ASSERT(size(RDM_2) == nAcpr2)
        n_lines = inv_triang_number(nAcpr2)

        file_id = arbitrary_magic_number
        file_id = isfreeunit(file_id)
        i = 1
        call molcas_open(file_id, trim(path))
            do curr_line = 1, n_lines
                read(file_id, *, iostat=io_err) RDM_2(i : i + curr_line - 1)
                call verify_(io_err == 0, 'Error on reading 2-RDMs.')
                i = i + curr_line
            end do
        close(file_id)
    end subroutine

    pure integer function triangular_number(n)
        integer, intent(in) :: n
        triangular_number = n * (n + 1) / 2
    end function

      !> This is the inverse function of triangular_number
    pure function inv_triang_number(n) result(res)
        integer, intent(in) :: n
        integer :: res
        res = nint(-0.5_wp + sqrt(0.25_wp + real(2 * n, kind=wp)))
    end function


    subroutine write_RDM(RDM, i_unit)
        real(wp), intent(in) :: RDM(:)
        integer, intent(in) :: i_unit

        integer :: io_err, curr_line, i, n_lines, j

        n_lines = inv_triang_number(size(RDM))

        i = 1
        do curr_line = 1, n_lines
            do j = i, i + curr_line - 1
                write(i_unit, '(E25.15)', advance='no', iostat=io_err) RDM(j)
                call verify_(io_err == 0, 'Error on writing RDM.')
            end do
            write(i_unit, *)
            i = i + curr_line
        end do
    end subroutine

end module CC_CI_mod
