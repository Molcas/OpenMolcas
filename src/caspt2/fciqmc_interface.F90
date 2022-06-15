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
      real(wp), intent(out) :: g1(nLev, nLev), g2(nLev, nLev, nLev, nLev), g3(*), &
                               f1(nLev, nLev), f2(nLev, nLev, nLev, nLev), f3(*)
      integer(1), intent(in) :: idxG3(6, *)
      logical :: proceed_found
      character(len=1024) :: WorkDir
      integer :: err

      proceed_found = .false.
      call getcwd_(WorkDir, err)
      write(u6, '(a)') 'Waiting for 3RDM and contracted Fock matrix.'
      write(u6, '(a)') 'Copy the file "fciqmc.caspt2.' // str(mstate(jState)) // &
        &'.h5" from M7 or NECI into the run directory:'
      write(u6, '(a)') trim(WorkDir)
      write(u6, '(a)') 'Afterwards, create a file "PROCEED" in the same folder.'
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
                 len6index(2), i, t, u, v, x, y, z
      logical :: tExist
      integer, allocatable :: indices(:,:)
      real(wp), allocatable :: values(:)
      real(wp) :: f3_temp(nLev, nLev, nLev, nLev, nLev, nLev), &
                  g3_temp(nLev, nLev, nLev, nLev, nLev, nLev)

      call f_Inquire('fciqmc.caspt2.' // str(iroot) // '.h5', tExist)
      call verify_(tExist, 'fciqmc.caspt2.' // str(iroot) // '.h5 does not exist.')
      hdf5_file = mh5_open_file_r('fciqmc.caspt2.' // str(iroot) // '.h5')
      hdf5_group = mh5_open_group(hdf5_file, 'archive/rdms/3300')
      hdf5_dset = mh5_open_dset(hdf5_group, 'indices')
      len6index(:) = 0
      call mh5_get_dset_dims(hdf5_dset, len6index)
      ! The HDF5 utilities transfer the indices and values as is,
      ! i.e. in C-style row major and with array indices -1 wrt. Fortran.
      ! When we take the HDF5 from NECI instead of M7 we have to take note.
      call mma_allocate(indices, 6, len6index(2))  ! six indices p, q, r, s, t, u
      call mma_allocate(values, len6index(2))
      indices(:,:) = 0
      values(:) = 0.0_wp
      call mh5_fetch_dset(hdf5_group, 'values', values)
      call mh5_fetch_dset(hdf5_group, 'indices', indices)
      call mh5_close_group(hdf5_group)

      g3_temp(:,:,:,:,:,:) = 0.0_wp
      do i = 1, len6index(2)
        ! refer to comment above
        t = indices(1,i) + 1; u = indices(2,i) + 1; v = indices(3,i) + 1
        x = indices(4,i) + 1; y = indices(5,i) + 1; z = indices(6,i) + 1
        g3_temp(t,u,v,x,y,z) = values(i)
        call apply_12fold_symmetry(g3_temp, t, u, v, x, y, z)
      end do
      call mma_deallocate(indices)
      call mma_deallocate(values)
      do i = 1, nG3
        t = idxG3(1,i); u = idxG3(2,i); v = idxG3(3,i)
        x = idxG3(4,i); y = idxG3(5,i); z = idxG3(6,i)
        g3(i) = g3_temp(t,u,v,x,y,z)
      end do
      write(u6,'(a)') "Completed the 3RDM transfer."

      hdf5_group = mh5_open_group(hdf5_file, 'archive/rdms/4400f')
      hdf5_dset = mh5_open_dset(hdf5_group, 'indices')
      len6index(:) = 0
      call mh5_get_dset_dims(hdf5_dset, len6index)
      ! refer to 3RDM
      call mma_allocate(indices, 6, len6index(2))  ! six indices p, q, r, s, t, u
      call mma_allocate(values, len6index(2))
      indices(:,:) = 0
      values(:) = 0.0_wp
      call mh5_fetch_dset(hdf5_group, 'values', values)
      call mh5_fetch_dset(hdf5_group, 'indices', indices)
      call mh5_close_group(hdf5_group)

      f3_temp(:,:,:,:,:,:) = 0.0_wp
      do i = 1, len6index(2)
        ! refer to 3RDM
        t = indices(1,i) + 1; u = indices(2,i) + 1; v = indices(3,i) + 1
        x = indices(4,i) + 1; y = indices(5,i) + 1; z = indices(6,i) + 1
        f3_temp(t,u,v,x,y,z) = values(i)
        call apply_12fold_symmetry(f3_temp, t, u, v, x, y, z)
      end do
      call mma_deallocate(indices)
      call mma_deallocate(values)
      do i = 1, nG3
        t = idxG3(1,i); u = idxG3(2,i); v = idxG3(3,i)
        x = idxG3(4,i); y = idxG3(5,i); z = idxG3(6,i)
        f3(i) = f3_temp(t,u,v,x,y,z)
      end do
      write(u6,'(a)') "Completed the F4RDM transfer."

      call calc_f2_and_g2(nActel, nLev, f3_temp, g3_temp, f2, g2)
      write(u6,'(a)') "Computed F2 and G2."
      call calc_f1_and_g1(nActel, nLev, f2, g2, f1, g1)
      write(u6,'(a)') "Computed F1 and G1."

      call mh5_close_file(hdf5_file)

      contains

        pure subroutine apply_12fold_symmetry(array, t, u, v, x, y, z)
          real(wp), intent(inout) :: array(:,:,:,:,:,:)
          integer, intent(in) :: t, u, v, x, y, z

          array(t,u,y,z,v,x) = array(t,u,v,x,y,z)
          array(u,t,x,v,z,y) = array(t,u,v,x,y,z)
          array(u,t,z,y,x,v) = array(t,u,v,x,y,z)
          array(v,x,t,u,y,z) = array(t,u,v,x,y,z)
          array(v,x,y,z,t,u) = array(t,u,v,x,y,z)
          array(x,v,u,t,z,y) = array(t,u,v,x,y,z)
          array(x,v,z,y,u,t) = array(t,u,v,x,y,z)
          array(y,z,t,u,v,x) = array(t,u,v,x,y,z)
          array(y,z,v,x,t,u) = array(t,u,v,x,y,z)
          array(z,y,u,t,x,v) = array(t,u,v,x,y,z)
          array(z,y,x,v,u,t) = array(t,u,v,x,y,z)
        end subroutine apply_12fold_symmetry

        pure subroutine calc_f2_and_g2(nAct, nLev, f3_temp, g3_temp, f2, g2)
          integer, intent(in) :: nAct, nLev
          real(wp), intent(in) :: f3_temp(nLev, nLev, nLev, nLev, nLev, nLev), &
                                  g3_temp(nLev, nLev, nLev, nLev, nLev, nLev)
          real(wp), intent(inout) :: f2(nLev, nLev, nLev, nLev), &
                                     g2(nLev, nLev, nLev, nLev)
          integer :: t, u, v, x, w

          f2(:,:,:,:) = 0.0_wp
          g2(:,:,:,:) = 0.0_wp
          do w = 1, nLev
            do x = 1, nLev
              do v = 1, nLev
                do u = 1, nLev
                  do t = 1, nLev
                    f2(t,u,v,x) = f2(t,u,v,x) + f3_temp(t,u,v,x,w,w)
                    g2(t,u,v,x) = g2(t,u,v,x) + g3_temp(t,u,v,x,w,w)
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
          integer :: t, u, w

          f1(:,:) = 0.0_wp
          g1(:,:) = 0.0_wp
          do w = 1, nLev
            do u = 1, nLev
              do t = 1, nLev
                f1(t,u) = f1(t,u) + f2(t,u,w,w)
                g1(t,u) = g1(t,u) + g2(t,u,w,w)
              end do
            end do
          end do
          f1(:,:) = f1(:,:) / (nAct - 2)
          g1(:,:) = g1(:,:) / (nAct - 1)
        end subroutine calc_f1_and_g1

    end subroutine load_fciqmc_mats

end module fciqmc_interface
