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
! Copyright (C) 2022-2023, Arta Safari                                 *
!***********************************************************************

module fciqmc_interface

#ifdef _MOLCAS_MPP_
use mpi, only: MPI_COMM_WORLD, MPI_LOGICAL
use Para_Info, only: Is_Real_Par
use Definitions, only: MPIInt
#endif
#ifdef _HDF5_
use mh5, only: mh5_close_file, mh5_close_group, mh5_fetch_dset, mh5_get_dset_dims, mh5_open_dset, mh5_open_file_r, mh5_open_group
#endif
use Para_Info, only: MyRank
use linalg_mod, only: verify_
use fortran_strings, only: str
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One
use, intrinsic :: iso_fortran_env, only: int8
use Definitions, only: wp, iwp, u6

implicit none
private

logical(kind=iwp) :: DoFCIQMC = .false., NonDiagonal = .false.

public :: DoFCIQMC, load_fciqmc_g1, mkfg3fciqmc, NonDiagonal

#include "macros.fh"

contains

!>  @brief
!>    Load 1RDM into poly1 in order to compute Fock matrix.
!>    Also serves as barrier for user to supply the other
!>    CASPT2 intermediates.
!>
!>  @param[in]     nLev      number of levels
!>  @param[inout]  g1        dense redundant 1RDM
!>  @param[in]     iroot     CASSCF root number
subroutine load_fciqmc_g1(nLev,g1,iroot)

  use filesystem, only: getcwd_
  use caspt2_global, only: mState, jState

  integer(kind=iwp), intent(in) :: nLev, iroot
  real(kind=wp), intent(inout) :: g1(nLev,nLev)
# ifdef _HDF5_
  integer(kind=iwp) :: err, hdf5_dset, hdf5_file, hdf5_group, i, len2index(2), t, u
  logical(kind=iwp) :: proceed_found, tExist
  character(len=1024) :: WorkDir
  integer(kind=iwp), allocatable :: indices(:,:)
  real(kind=wp), allocatable :: values(:)
# ifdef _MOLCAS_MPP_
  integer(kind=MPIInt) :: error
  integer(kind=MPIInt), parameter :: ROOT = 0_MPIInt
# endif

  proceed_found = .false.
  call getcwd_(WorkDir,err)
  write(u6,'(4x,a)') 'Waiting for the 3RDM and contracted Fock matrix.'
  write(u6,'(4x,a)') 'First copy the required files into the M7 work directory:'
  if (NonDiagonal) then
    write(u6,'(8x,a)') 'cp '//trim(WorkDir)//'/fockdump.h5 $M7_WORKDIR'
    write(u6,'(4x,a)') 'Use the same FciDump as for the preceeding CASCI.'
  else
    write(u6,'(8x,a)') 'cp '//trim(WorkDir)//'/{fockdump.h5,caspt2.FciDmp.h5} $M7_WORKDIR'
  end if
  write(u6,'(4x,a)') 'With these files run the FCIQMC dynamic.'
  write(u6,'(4x,a)') 'Copy the file M7.rdm.h5 as "fciqmc.caspt2.'//str(mstate(jState))//'.h5" into the run directory.'
  write(u6,'(4x,a)') 'Afterwards, create a file "PROCEED" in the same folder:'
  write(u6,'(8x,a)') 'cp $M7_WORKDIR/M7.rdm.h5 '//trim(WorkDir)//'/fciqmc.caspt2.'//str(mstate(jState))//'.h5'
  write(u6,'(8x,a)') 'touch '//trim(WorkDir)//'/PROCEED'

  do while (.not. proceed_found)
    call sleepf(1)
    if (myrank == 0) call f_Inquire('PROCEED',proceed_found)
#   ifdef _MOLCAS_MPP_
    if (is_real_par()) call MPI_Bcast(proceed_found,1_MPIInt,MPI_LOGICAL,ROOT,MPI_COMM_WORLD,error)
#   endif
  end do
  if (myRank == 0) then
    write(u6,'(a)') 'PROCEED file found. Continuing with CASPT2.'
  else
    call bcast_2RDM('fciqmc.caspt2.'//str(iroot)//'.h5')
  end if

  call f_Inquire('fciqmc.caspt2.'//str(iroot)//'.h5',tExist)
  call verify_(tExist,'fciqmc.caspt2.'//str(iroot)//'.h5 does not exist.')
  hdf5_file = mh5_open_file_r('fciqmc.caspt2.'//str(iroot)//'.h5')
  hdf5_group = mh5_open_group(hdf5_file,'/spinfree/1100')
  hdf5_dset = mh5_open_dset(hdf5_group,'indices')
  len2index(:) = 0
  call mh5_get_dset_dims(hdf5_dset,len2index)
  call mma_allocate(indices,2,len2index(2))
  call mma_allocate(values,len2index(2))
  indices(:,:) = 0
  values(:) = Zero
  call mh5_fetch_dset(hdf5_group,'values',values)
  call mh5_fetch_dset(hdf5_group,'indices',indices)
  call mh5_close_group(hdf5_group)
  g1(:,:) = Zero
  do i=1,len2index(2)
    t = indices(1,i)+1
    u = indices(2,i)+1
    g1(t,u) = values(i)
  end do
  if (NonDiagonal) then
    call transform_1rdm(g1,nLev)
    write(u6,'(a)') 'Transformed 1RDM to pseudo-canonical orbitals.'
  end if
  write(u6,'(a)') 'Completed the 1RDM transfer.'
  call mma_deallocate(indices)
  call mma_deallocate(values)
  call mh5_close_file(hdf5_file)

contains

  !>  @brief
  !>    Transform 1RDM to pseudo-canonical orbitals. To this end,
  !>    read Fock matrix eigenvectors from fockdump.h5.
  !>
  !>  @param[inout]    g1        dense redundant 1RDM
  !>  @param[in]       nLev      number of levels
  subroutine transform_1rdm(g1,nLev)

    integer(kind=iwp), intent(in) :: nLev
    real(kind=wp), intent(inout) :: g1(nLev,nLev)
    integer(kind=iwp) :: hdf5_file, hdf5_group
    real(kind=wp) :: fockvecs(nLev,nLev)
    logical(kind=iwp) :: tExist

    if (myRank /= 0) call bcast_2RDM('fockdump.h5')
    call f_Inquire('fockdump.h5',tExist)
    call verify_(tExist,'fockdump.h5 does not exist.')
    hdf5_file = mh5_open_file_r('fockdump.h5')
    hdf5_group = mh5_open_group(hdf5_file,'/')
    call mh5_fetch_dset(hdf5_group,'ACT_FOCK_EIGVECS',fockvecs)
    call mh5_close_group(hdf5_group)
    call transmat(g1,fockvecs,nLev)
    call mh5_close_file(hdf5_file)

  end subroutine transform_1rdm

# else
  unused_var(nLev)
  unused_var(g1)
  unused_var(iroot)
# endif

end subroutine load_fciqmc_g1

!>  @brief
!>    Wrapper to collect all required density and Fock matrices and feed
!>    them into poly3. Interface consistent with caspt2 mkfg3.f
!>
!>  @param[out]    g1        dense redundant 1RDM
!>  @param[out]    g2        dense redundant 2RDM
!>  @param[out]    g3        sparse 3RDM
!>  @param[out]    f1        dense contraction of Fockian with 2RDM
!>  @param[out]    f2        dense contraction of Fockian with 3RDM
!>  @param[out]    f3        sparse contraction of Fockian with 4RDM
!>  @param[in]     idxG3     Table containing the active space indices
subroutine mkfg3fciqmc(g1,g2,g3,f1,f2,f3,idxG3,nLev)

# ifdef _HDF5_
  use caspt2_global, only: jState, mState, nG3
# endif

  integer(kind=iwp), intent(in) :: nLev
  real(kind=wp), intent(inout) :: g1(nLev,nLev), g2(nLev,nLev,nLev,nLev), g3(*), f1(nLev,nLev), f2(nLev,nLev,nLev,nLev), f3(*)
  integer(kind=int8), intent(in) :: idxG3(6,*)

# ifdef _HDF5_
  call load_fciqmc_mats(nLev,idxG3,nG3,g3,g2,g1,f3,f2,f1,mstate(jState))
# else
  unused_var(idxG3(1,1))
  unused_var(g3(1))
  unused_var(f3(1))
  unused_var(g2)
  unused_var(g1)
  unused_var(f2)
  unused_var(f1)
# endif

end subroutine mkfg3fciqmc

#ifdef _HDF5_
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
subroutine load_fciqmc_mats(nLev,idxG3,nG3,g3,g2,g1,f3,f2,f1,iroot)

  use caspt2_global, only: nActEl

  integer(kind=iwp), intent(in) :: nLev, nG3, iroot
  integer(kind=int8), intent(in) :: idxG3(6,nG3)
  real(kind=wp), intent(inout) :: g3(*), g2(nLev,nLev,nLev,nLev), g1(nLev,nLev), f3(*), f2(nLev,nLev,nLev,nLev), f1(nLev,nLev)
  integer(kind=iwp) :: hdf5_dset, hdf5_file, hdf5_group, i, len6index(2), t, u, v, x, y, z
  real(kind=wp) :: f3_temp(nLev,nLev,nLev,nLev,nLev,nLev), g3_temp(nLev,nLev,nLev,nLev,nLev,nLev)
  logical(kind=iwp) :: tExist
  integer(kind=iwp), allocatable :: indices(:,:)
  real(kind=wp), allocatable :: values(:)
# ifdef _DEBUGPRINT_
  real(kind=wp) :: cpu, cpu0, cpu1, tio, tio0, tio1, trace
# endif

  if (myRank /= 0) call bcast_2RDM('fciqmc.caspt2.'//str(iroot)//'.h5')
  call f_Inquire('fciqmc.caspt2.'//str(iroot)//'.h5',tExist)
  call verify_(tExist,'fciqmc.caspt2.'//str(iroot)//'.h5 does not exist.')
  hdf5_file = mh5_open_file_r('fciqmc.caspt2.'//str(iroot)//'.h5')
  hdf5_group = mh5_open_group(hdf5_file,'/spinfree/3300')
  hdf5_dset = mh5_open_dset(hdf5_group,'indices')
  len6index(:) = 0
  call mh5_get_dset_dims(hdf5_dset,len6index)
  call mma_allocate(indices,6,len6index(2))
  call mma_allocate(values,len6index(2))
  indices(:,:) = 0
  values(:) = Zero
  call mh5_fetch_dset(hdf5_group,'values',values)
  call mh5_fetch_dset(hdf5_group,'indices',indices)
  call mh5_close_group(hdf5_group)
  g3_temp(:,:,:,:,:,:) = Zero
  do i=1,len6index(2)
    t = indices(1,i)+1
    u = indices(4,i)+1
    v = indices(2,i)+1
    x = indices(5,i)+1
    y = indices(3,i)+1
    z = indices(6,i)+1
    call apply_12fold_symmetry(g3_temp,t,u,v,x,y,z,values(i))
  end do
  call mma_deallocate(indices)
  call mma_deallocate(values)

# ifdef _DEBUGPRINT_
  trace = Zero
  do v=1,nLev
    do u=1,nLev
      do t=1,nLev
        trace = trace+g3_temp(t,t,u,u,v,v)
      end do
    end do
  end do
  write(u6,'(a,f12.5)') 'Trace 3RDM: ',trace
# endif

  if (NonDiagonal) then
#   ifdef _DEBUGPRINT_
    call timing(cpu0,cpu,tio0,tio)
#   endif
    write(u6,'(a)') 'Transformed 3RDM to pseudo-canonical orbitals.'
    call transform_six_index(g3_temp,nLev)
#   ifdef _DEBUGPRINT_
    call timing(cpu1,cpu,tio1,tio)
    write(u6,*) 'Wall time 3RDM transform: ',tio1-tio0
    trace = Zero
    do v=1,nLev
      do u=1,nLev
        do t=1,nLev
          trace = trace+g3_temp(t,t,u,u,v,v)
        end do
      end do
    end do
    write(u6,'(a,f12.5)') 'Trace transformed 3RDM: ',trace
#   endif
  end if

  do i=1,nG3
    t = idxG3(1,i)
    u = idxG3(2,i)
    v = idxG3(3,i)
    x = idxG3(4,i)
    y = idxG3(5,i)
    z = idxG3(6,i)
    g3(i) = g3_temp(t,u,v,x,y,z)
  end do
  write(u6,'(a)') 'Completed the 3RDM transfer.'

  hdf5_group = mh5_open_group(hdf5_file,'/spinfree/4400f')
  hdf5_dset = mh5_open_dset(hdf5_group,'indices')
  len6index(:) = 0
  call mh5_get_dset_dims(hdf5_dset,len6index)
  call mma_allocate(indices,6,len6index(2))
  call mma_allocate(values,len6index(2))
  indices(:,:) = 0
  values(:) = Zero
  call mh5_fetch_dset(hdf5_group,'values',values)
  call mh5_fetch_dset(hdf5_group,'indices',indices)
  call mh5_close_group(hdf5_group)
  f3_temp(:,:,:,:,:,:) = Zero
  do i=1,len6index(2)
    t = indices(1,i)+1
    u = indices(4,i)+1
    v = indices(2,i)+1
    x = indices(5,i)+1
    y = indices(3,i)+1
    z = indices(6,i)+1
    call apply_12fold_symmetry(f3_temp,t,u,v,x,y,z,values(i))
  end do
  call mma_deallocate(indices)
  call mma_deallocate(values)

  if (NonDiagonal) then
    call transform_six_index(f3_temp,nLev)
    write(u6,'(a)') 'Transformed F.4RDM to pseudo-canonical orbitals.'
  end if

  do i=1,nG3
    t = idxG3(1,i)
    u = idxG3(2,i)
    v = idxG3(3,i)
    x = idxG3(4,i)
    y = idxG3(5,i)
    z = idxG3(6,i)
    f3(i) = f3_temp(t,u,v,x,y,z)
  end do
  write(u6,'(a)') 'Completed the F.4RDM transfer.'

  call calc_f2_and_g2(nActel,nLev,f3_temp,g3_temp,f2,g2)
  write(u6,'(a)') 'Computed F2 and G2.'
  call calc_f1_and_g1(nActel,nLev,f2,g2,f1,g1)
  write(u6,'(a)') 'Computed F1 and G1.'

  call mh5_close_file(hdf5_file)

contains

  !>  @brief
  !>    Transform 3RDM and F.4RDM to pseudo-canonical orbitals.
  !>
  !>  @param[inout]    g3        dense redundant 3RDM
  !>  @param[inout]    f3        dense redundant F.4RDM
  !>  @param[in]       nLev      number of levels
  subroutine transform_six_index(six_index,nLev)

    integer(kind=iwp), intent(in) :: nLev
    real(kind=wp), intent(inout) :: six_index(nLev,nLev,nLev,nLev,nLev,nLev)
    integer(kind=iwp) :: hdf5_file, hdf5_group, iter, u2, v2, x2, y2, z2  ! prevent sharing scope with upper function
    logical(kind=iwp) :: tExist
    real(kind=wp) :: buffer(nLev), buffer2(nLev), fockvecs(nLev,nLev)

    call f_Inquire('fockdump.h5',tExist)
    call verify_(tExist,'fockdump.h5 does not exist.')
    hdf5_file = mh5_open_file_r('fockdump.h5')
    hdf5_group = mh5_open_group(hdf5_file,'/')
    call mh5_fetch_dset(hdf5_group,'ACT_FOCK_EIGVECS',fockvecs)
    call mh5_close_group(hdf5_group)

    buffer(:) = Zero
    buffer2(:) = Zero
    do iter=1,6
      do z2=1,nLev
        do y2=1,nLev
          do x2=1,nLev
            do v2=1,nLev
              do u2=1,nLev
                select case (iter)
                  case (1)
                    buffer(:) = six_index(:,u2,v2,x2,y2,z2)
                    call dgemv_('T',nLev,nLev,One,fockvecs,nLev,buffer,1,Zero,buffer2,1)
                    six_index(:,u2,v2,x2,y2,z2) = buffer2(:)
                  case (2)
                    buffer(:) = six_index(u2,:,v2,x2,y2,z2)
                    call dgemv_('T',nLev,nLev,One,fockvecs,nLev,buffer,1,Zero,buffer2,1)
                    six_index(u2,:,v2,x2,y2,z2) = buffer2(:)
                  case (3)
                    buffer(:) = six_index(u2,v2,:,x2,y2,z2)
                    call dgemv_('T',nLev,nLev,One,fockvecs,nLev,buffer,1,Zero,buffer2,1)
                    six_index(u2,v2,:,x2,y2,z2) = buffer2(:)
                  case (4)
                    buffer(:) = six_index(u2,v2,x2,:,y2,z2)
                    call dgemv_('T',nLev,nLev,One,fockvecs,nLev,buffer,1,Zero,buffer2,1)
                    six_index(u2,v2,x2,:,y2,z2) = buffer2(:)
                  case (5)
                    buffer(:) = six_index(u2,v2,x2,y2,:,z2)
                    call dgemv_('T',nLev,nLev,One,fockvecs,nLev,buffer,1,Zero,buffer2,1)
                    six_index(u2,v2,x2,y2,:,z2) = buffer2(:)
                  case (6)
                    buffer(:) = six_index(u2,v2,x2,y2,z2,:)
                    call dgemv_('T',nLev,nLev,One,fockvecs,nLev,buffer,1,Zero,buffer2,1)
                    six_index(u2,v2,x2,y2,z2,:) = buffer2(:)
                end select
              end do
            end do
          end do
        end do
      end do
    end do
    call mh5_close_file(hdf5_file)

  end subroutine transform_six_index

  pure subroutine apply_12fold_symmetry(array,t,u,v,x,y,z,val)
    ! G3 has 12 permutational symmetries, since the spin indices of
    ! the (t,u), (v,x) and (y,z) indices have to match up.

    real(kind=wp), intent(inout) :: array(:,:,:,:,:,:)
    integer(kind=iwp), intent(in) :: t, u, v, x, y, z
    real(kind=wp), intent(in) :: val

    array(t,u,v,x,y,z) = val
    array(t,u,y,z,v,x) = val
    array(v,x,t,u,y,z) = val
    array(v,x,y,z,t,u) = val
    array(y,z,t,u,v,x) = val
    array(y,z,v,x,t,u) = val
    array(u,t,x,v,z,y) = val
    array(u,t,z,y,x,v) = val
    array(x,v,u,t,z,y) = val
    array(x,v,z,y,u,t) = val
    array(z,y,u,t,x,v) = val
    array(z,y,x,v,u,t) = val

  end subroutine apply_12fold_symmetry

  pure subroutine calc_f2_and_g2(nAct,nLev,f3_temp,g3_temp,f2,g2)

    integer(kind=iwp), intent(in) :: nAct, nLev
    real(kind=wp), intent(in) :: f3_temp(nLev,nLev,nLev,nLev,nLev,nLev), g3_temp(nLev,nLev,nLev,nLev,nLev,nLev)
    real(kind=wp), intent(inout) :: f2(nLev,nLev,nLev,nLev), g2(nLev,nLev,nLev,nLev)
    integer(kind=iwp) :: w

    f2(:,:,:,:) = Zero
    g2(:,:,:,:) = Zero
    do w=1,nLev
      f2(:,:,:,:) = f2(:,:,:,:)+f3_temp(:,:,:,:,w,w)
      g2(:,:,:,:) = g2(:,:,:,:)+g3_temp(:,:,:,:,w,w)
    end do
    f2(:,:,:,:) = f2(:,:,:,:)/(nAct-3)
    g2(:,:,:,:) = g2(:,:,:,:)/(nAct-2)

  end subroutine calc_f2_and_g2

  pure subroutine calc_f1_and_g1(nAct,nLev,f2,g2,f1,g1)

    integer(kind=iwp), intent(in) :: nAct, nLev
    real(kind=wp), intent(in) :: f2(nLev,nLev,nLev,nLev), g2(nLev,nLev,nLev,nLev)
    real(kind=wp), intent(inout) :: f1(nLev,nLev), g1(nLev,nLev)
    integer(kind=iwp) :: w

    f1(:,:) = Zero
    g1(:,:) = Zero
    do w=1,nLev
      f1(:,:) = f1(:,:)+f2(:,:,w,w)
      g1(:,:) = g1(:,:)+g2(:,:,w,w)
    end do
    f1(:,:) = f1(:,:)/(nAct-2)
    g1(:,:) = g1(:,:)/(nAct-1)

  end subroutine calc_f1_and_g1

end subroutine load_fciqmc_mats

subroutine bcast_2RDM(InFile)

  use filesystem, only: get_errno_, strerror_, symlink_

  character(len=*), intent(in) :: InFile
  character(len=1024) :: master
  integer(kind=iwp) :: err, lmaster1

  call prgmtranslate_master(InFile,master,lmaster1)
  call symlink_(trim(master),trim(InFile),err)
  if (err == 0) write(u6,*) strerror_(get_errno_())

end subroutine bcast_2RDM

#endif

end module fciqmc_interface
