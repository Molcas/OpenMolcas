!**********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!                                                                      *
! Copyright (C) 2022, Arta Safari                                      *
!***********************************************************************
module fciqmc_interface

#ifdef _MOLCAS_MPP_
    use mpi
    use Para_Info, only: Is_Real_Par
    use definitions, only: MPIInt
#endif
    use Para_Info, only: MyRank
    use filesystem, only : getcwd_
    use definitions, only: wp, u6
    use linalg_mod, only: verify_, abort_
    use fortran_strings, only: str
    use stdalloc, only: mma_allocate, mma_deallocate
    use mh5, only: mh5_open_file_r, mh5_close_file, &
                   mh5_open_group, mh5_close_group, &
                   mh5_open_dset, mh5_close_dset, mh5_fetch_dset, mh5_get_dset_dims, &
                   mh5_exists_dset

    implicit none

    private
    public :: mkfg3fciqmc, DoFCIQMC
    logical :: DoFCIQMC = .false.

    contains

    !>  @brief
    !>    Wrapper to collect all required density and Fock matrices and feed
    !>    them into poly3.
    !>
    !>  @param[out]    g1        dense redundant 1RDM
    !>  @param[out]    g2        dense redundant 2RDM
    !>  @param[out]    g3        sparse 3RDM
    !>  @param[out]    f1        dense contraction of Fockian with 2RDM
    !>  @param[out]    f2        dense contraction of Fockian with 3RDM
    !>  @param[out]    f3        sparse contraction of Fockian with 4RDM
    !>  @param[in]     idxG3     Table containing the active space indices
    subroutine mkfg3fciqmc(g1, g2, g3, f1, f2, f3, idxG3)
#ifdef NAGFOR
        use f90_unix_proc, only: sleep
#endif
        use caspt2_data, only: nG3, nLev, mState, jState
        real(wp), intent(inout) :: g1(nLev, nLev), g2(nLev, nLev, nLev, nLev), g3(*), &
                                 f1(nLev, nLev), f2(nLev, nLev, nLev, nLev), f3(*)
        integer(1), intent(in) :: idxG3(6, *)
        logical :: proceed_found
        character(len=1024) :: WorkDir
        integer :: err

        proceed_found = .false.
        call getcwd_(WorkDir, err)
        write(u6, '(4x,a)') 'Waiting for the 3RDM and contracted Fock matrix.'
        write(u6, '(4x,a)') 'Copy the file "fciqmc.caspt2.' // str(mstate(jState)) // &
            &'.h5" from M7 into the run directory:'
        write(u6, '(4x,a)') 'Afterwards, create a file "PROCEED" in the same folder.'
        write(u6, '(8x,a)') 'cp fciqmc.caspt2.'// str(mstate(jState)) //'.h5 ' // trim(WorkDir)
        write(u6, '(8x,a)') 'touch ' // trim(WorkDir) // '/PROCEED'
        do while(.not. proceed_found)
            call sleep(1)
            if (myrank == 0) call f_Inquire('PROCEED', proceed_found)
#ifdef _MOLCAS_MPP_
            if (is_real_par()) then
                call MPI_Bcast(proceed_found, 1_MPIInt, MPI_LOGICAL, &
                               ROOT, MPI_COMM_WORLD, error)
            end if
#endif
        end do
        if (myrank == 0) then
            write(u6, '(a)') 'PROCEED file found. Continuing with CASPT2.'
        end if

        call load_fciqmc_mats(nLev, idxG3, nG3, g3, g2, g1, &
                                f3, f2, f1, mstate(jState))
    end subroutine mkfg3fciqmc


    !>  @brief
    !>    Read stochastically sampled 3RDMs and contracted Fock tensors
    !>    stored in HDF5 format.
    !>
    !>  @param[in]     nLev       Number of Levels in the GUGA formalism
    !>  @param[in]     idxG3      Table containing the active space indices
    !>  @param[in]     nG3        Number of 3RDM elements
    !>  @param[in]     g3         3RDM
    !>  @param[in]     g2         2RDM
    !>  @param[in]     g1         1RDM
    !>  @param[in]     f3         contracted Fock matrix with 4RDM
    !>  @param[in]     f2         contracted Fock matrix with 3RDM
    !>  @param[in]     f1         contracted Fock matrix with 2RDM
    !>  @param[in]     iroot      MCSCF root number.
    subroutine load_fciqmc_mats(nLev, idxG3, nG3, g3, g2, g1, f3, f2, f1, iroot)
        use caspt2_data, only: nActEl
        integer, intent(in) :: nLev
        integer(1), intent(in) :: idxG3(6, nG3)
        integer, intent(in) :: nG3
        real(wp), intent(inout) :: g3(*), g2(nLev, nLev, nLev, nLev), g1(nLev, nLev), &
                                   f3(*), f2(nLev, nLev, nLev, nLev), f1(nLev, nLev)
        integer, intent(in) :: iroot
        integer :: hdf5_file, hdf5_group, hdf5_dset, &
                   len6index(2), i, p, q, r, s, t, u
        logical :: tExist
        integer, allocatable :: indices(:,:)
        real(wp), allocatable :: values(:)
        real(wp) :: f3_temp(nLev, nLev, nLev, nLev, nLev, nLev), &
                    g3_temp(nLev, nLev, nLev, nLev, nLev, nLev)

        call f_Inquire('fciqmc.caspt2.' // str(iroot) // '.h5', tExist)
        call verify_(tExist, 'fciqmc.caspt2.' // str(iroot) // '.h5 does not exist.')
        hdf5_file = mh5_open_file_r('fciqmc.caspt2.' // str(iroot) // '.h5')
        hdf5_group = mh5_open_group(hdf5_file, 'archive/rdms/sf_3300')
        hdf5_dset = mh5_open_dset(hdf5_group, 'indices')
        len6index(:) = 0
        call mh5_get_dset_dims(hdf5_dset, len6index)
        call mma_allocate(indices, 6, len6index(2))
        call mma_allocate(values, len6index(2))
        indices(:,:) = 0
        values(:) = 0.0_wp
        call mh5_fetch_dset(hdf5_group, 'values', values)
        call mh5_fetch_dset(hdf5_group, 'indices', indices)
        call mh5_close_group(hdf5_group)
        g3_temp(:,:,:,:,:,:) = 0.0_wp
        do i = 1, len6index(2)
            ! The HDF5 utilities transfer the indices and values as is,
            ! i.e. in C-style row major and with array indices -1 wrt. Fortran.
            ! When we take the HDF5 from NECI instead of M7 we have to remove the +1.
            p = indices(1,i) + 1; q = indices(2,i) + 1; r = indices(3,i) + 1
            s = indices(4,i) + 1; t = indices(5,i) + 1; u = indices(6,i) + 1
            ! note the index ordering
            call apply_12fold_symmetry(g3_temp, p, q, r, s, t, u, values(i))
        end do
        call mma_deallocate(indices)
        call mma_deallocate(values)
        do i = 1, nG3
            p = idxG3(1,i); q = idxG3(2,i); r = idxG3(3,i)
            s = idxG3(4,i); t = idxG3(5,i); u = idxG3(6,i)
            g3(i) = g3_temp(p, r, t, q, s, u)
        end do
        write(u6,'(a)') "Completed the 3RDM transfer."

        hdf5_group = mh5_open_group(hdf5_file, 'archive/rdms/sf_4400f')
        hdf5_dset = mh5_open_dset(hdf5_group, 'indices')
        len6index(:) = 0
        call mh5_get_dset_dims(hdf5_dset, len6index)
        call mma_allocate(indices, 6, len6index(2))
        call mma_allocate(values, len6index(2))
        indices(:,:) = 0
        values(:) = 0.0_wp
        call mh5_fetch_dset(hdf5_group, 'values', values)
        call mh5_fetch_dset(hdf5_group, 'indices', indices)
        call mh5_close_group(hdf5_group)
        f3_temp(:,:,:,:,:,:) = 0.0_wp
        do i = 1, len6index(2)
            p = indices(1,i) + 1; q = indices(2,i) + 1; r = indices(3,i) + 1
            s = indices(4,i) + 1; t = indices(5,i) + 1; u = indices(6,i) + 1
            call apply_12fold_symmetry(f3_temp, p, q, r, s, t, u, values(i))
        end do
        call mma_deallocate(indices)
        call mma_deallocate(values)
        do i = 1, nG3
            p = idxG3(1,i); q = idxG3(2,i); r = idxG3(3,i)
            s = idxG3(4,i); t = idxG3(5,i); u = idxG3(6,i)
            f3(i) = f3_temp(p, r, t, q, s, u)
        end do
        write(u6,'(a)') "Completed the F4RDM transfer."

        call calc_f2_and_g2(nActel, nLev, f3_temp, g3_temp, f2, g2)
        write(u6,'(a)') "Computed F2 and G2."
        call calc_f1_and_g1(nActel, nLev, f2, g2, f1, g1)
        write(u6,'(a)') "Computed F1 and G1."

        call mh5_close_file(hdf5_file)

        contains

            pure subroutine apply_12fold_symmetry(array, p, q, r, s, t, u, val)
                ! These RDMs and Fock contractions all use the normal ordered definition
                ! and not the product-of-single-excitations one, i.e.:
                ! G3(p,q,r,s,t,u) = sum_{o1, o2, o3} < p_o1+ r_o2+ t_o3+ u_o3 s_o2 q_o1 >
                ! In the Gamma notation of doi: 10.1063/1.5140086 Appendix A this correspond to
                ! Gamma_prt,qsu.
                ! From this definition follows that G3 has 12 permutational symmetries, since
                ! the spin indices of the (p,q) (r,s) and (t,u) indices have to match up.
                real(wp), intent(inout) :: array(:,:,:,:,:,:)
                integer, intent(in) :: p, q, r, s, t, u
                real(wp), intent(in) :: val

                array(p, q, r, s, t, u) = val
                array(p, r, q, s, u, t) = val
                array(q, p, r, t, s, u) = val
                array(r, p, q, u, s, t) = val
                array(q, r, p, t, u, s) = val
                array(r, q, p, u, t, s) = val
                ! transpose the above
                array(s, t, u, p, q, r) = val
                array(s, u, t, p, r, q) = val
                array(t, s, u, q, p, r) = val
                array(u, s, t, r, p, q) = val
                array(t, u, s, q, r, p) = val
                array(u, t, s, r, q, p) = val
            end subroutine apply_12fold_symmetry

            pure subroutine calc_f2_and_g2(nAct, nLev, f3_temp, g3_temp, f2, g2)
                integer, intent(in) :: nAct, nLev
                real(wp), intent(in) :: f3_temp(nLev, nLev, nLev, nLev, nLev, nLev), &
                                        g3_temp(nLev, nLev, nLev, nLev, nLev, nLev)
                real(wp), intent(inout) :: f2(nLev, nLev, nLev, nLev), &
                                           g2(nLev, nLev, nLev, nLev)
                integer :: p, q, r, s, w

                f2(:,:,:,:) = 0.0_wp
                g2(:,:,:,:) = 0.0_wp
                do w = 1, nLev
                    do s = 1, nLev
                        do r = 1, nLev
                            do q = 1, nLev
                                do p = 1, nLev
                                    f2(p,q,r,s) = f2(p,q,r,s) + f3_temp(p,q,w,r,s,w)
                                    g2(p,q,r,s) = g2(p,q,r,s) + g3_temp(p,q,w,r,s,w)
                                end do
                            end do
                        end do
                    end do
                end do
                f2(:,:,:,:) = f2(:,:,:,:) / (nAct - 3)
                g2(:,:,:,:) = g2(:,:,:,:) / (nAct - 2)
            end subroutine calc_f2_and_g2

            pure subroutine calc_f1_and_g1(nAct, nLev, f2, g2, f1, g1)
                integer, intent(in) :: nAct, nLev
                real(wp), intent(in) :: f2(nLev, nLev, nLev, nLev), g2(nLev, nLev, nLev, nLev)
                real(wp), intent(inout) :: f1(nLev, nLev), g1(nLev, nLev)
                integer :: p, q, w

                f1(:,:) = 0.0_wp
                g1(:,:) = 0.0_wp
                do w = 1, nLev
                    do q = 1, nLev
                        do p = 1, nLev
                            ! this trace also has to switch indices, see above!
                            f1(p,q) = f1(p,q) + f2(p,w,q,w)
                            g1(p,q) = g1(p,q) + g2(p,w,q,w)
                        end do
                    end do
                end do
                f1(:,:) = f1(:,:) / (nAct - 2)
                g1(:,:) = g1(:,:) / (nAct - 1)
            end subroutine calc_f1_and_g1

    end subroutine load_fciqmc_mats

end module fciqmc_interface
