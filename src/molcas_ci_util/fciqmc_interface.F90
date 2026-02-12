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
! Copyright (C) 2022-2024, Arta Safari                                 *
!***********************************************************************

module fciqmc_interface

#ifdef _MOLCAS_MPP_
use mpi, only: MPI_COMM_WORLD, MPI_LOGICAL
use Para_Info, only: Is_Real_Par
use Definitions, only: MPIInt
#endif
use Definitions, only: wp, iwp, byte
#ifdef _HDF5_
use mh5, only: mh5_close_file, mh5_close_group, mh5_fetch_dset, mh5_get_dset_dims, mh5_open_dset, mh5_open_file_r, mh5_open_group
use Para_Info, only: MyRank
use caspt2_module, only: jstate, mstate, nActel
use pt2_guga, only: nG3
use linalg_mod, only: verify_
use fortran_strings, only: str
use filesystem, only: getcwd_
use Constants, only: Zero, One
use stdalloc, only: mma_allocate, mma_deallocate
use Definitions, only: u6
#endif

implicit none
private

logical(kind=iwp) :: DoFCIQMC = .false., NonDiagonal = .false., TransformToNormalOrder = .false.

public :: DoFCIQMC, load_fciqmc_g1, mkfg3fciqmc, NonDiagonal, TransformToNormalOrder

#include "macros.fh"

! new structure of GUGX module makes it necessary to pass nLev as function parameter every time :(

contains

!> @brief
!>   Load 1RDM into poly1 to compute the Fock matrix.
!>
!> @param[in]     nLev   number of levels
!> @param[inout]  g1     dense redundant 1RDM
!> @param[in]     iroot  CASSCF root number
subroutine load_fciqmc_g1(g1,iroot,nLev)

  integer(kind=iwp), intent(in) :: iroot, nLev
  real(kind=wp), intent(inout) :: g1(nLev,nLev)

# ifndef _HDF5_
  unused_var(nLev)
  unused_var(g1)
  unused_var(iroot)
# else
  call user_barrier()  ! copy required files into WorkDir
  call load_1RDM(g1)
  if (NonDiagonal) call transform_1rdm(g1)

contains

  !> @brief
  !>   Transform 1RDM to pseudo-canonical orbitals.
  subroutine transform_1rdm(g1)

    real(kind=wp), intent(inout) :: g1(nLev,nLev)
    real(kind=wp), allocatable :: fock_eigvecs(:,:), fockmat(:,:)

    call mma_allocate(fockmat,nLev,nLev,Label='fockmat')
    call mma_allocate(fock_eigvecs,nLev,nLev,Label='fock_eigvecs')
    call load_fockmat(fockmat,fock_eigvecs,nLev)
    call transmat(g1,fock_eigvecs,nLev)
    call mma_deallocate(fockmat)
    call mma_deallocate(fock_eigvecs)

  end subroutine transform_1rdm

  subroutine load_1RDM(g1)

    real(kind=wp), intent(inout) :: g1(nLev,nLev)
    integer(kind=iwp) :: hdf5_dset, hdf5_file, hdf5_group, i, len2index(2), t, u
    logical(kind=iwp) :: tExist
    integer(kind=iwp), allocatable :: indices(:,:)
    real(kind=wp), allocatable :: values(:)

    call f_Inquire('fciqmc.caspt2.'//str(iroot)//'.h5',tExist)
    call verify_(tExist,'fciqmc.caspt2.'//str(iroot)//'.h5 does not exist.')
    hdf5_file = mh5_open_file_r('fciqmc.caspt2.'//str(iroot)//'.h5')
    hdf5_group = mh5_open_group(hdf5_file,'/spinfree/1100/')
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
      g1(u,t) = values(i)
    end do
    call mma_deallocate(indices)
    call mma_deallocate(values)
    call mh5_close_file(hdf5_file)

  end subroutine load_1RDM

# endif

end subroutine load_fciqmc_g1

!> @brief
!>   Wrapper to collect all required density and Fock matrices and feed
!>   them into poly3. Interface consistent with caspt2 mkfg3.f
!>
!> @param[out]  g1     dense redundant 1RDM
!> @param[out]  g2     dense redundant 2RDM
!> @param[out]  g3     sparse 3RDM
!> @param[out]  f1     dense contraction of Fockian with 2RDM
!> @param[out]  f2     dense contraction of Fockian with 3RDM
!> @param[out]  f3     sparse contraction of Fockian with 4RDM
!> @param[in]   idxG3  Table containing the active space indices
subroutine mkfg3fciqmc(g1,g2,g3,f1,f2,f3,idxG3,nLev)

  integer(kind=iwp), intent(in) :: nLev
  real(kind=wp), intent(inout) :: g1(nLev,nLev), g2(nLev,nLev,nLev,nLev), g3(*), f1(nLev,nLev), f2(nLev,nLev,nLev,nLev), f3(*)
  integer(kind=byte), intent(in) :: idxG3(6,*)

# ifndef _HDF5_
  unused_var(g1)
  unused_var(g2)
  unused_var(g3(1))
  unused_var(f1)
  unused_var(f2)
  unused_var(f3(1))
  unused_var(idxG3(1,1))
# else
  call load_fciqmc_mats(idxG3,g3,g2,g1,f3,f2,f1,mstate(jState),nLev)
# endif

end subroutine mkfg3fciqmc

!> @brief
!>   Read stochastically sampled 3RDMs and contracted Fock tensors
!>   stored in HDF5 format.
!>
!> @param[in]  nLev   Number of Levels in the GUGA formalism
!> @param[in]  idxG3  Table containing the active space indices
!> @param[in]  nG3    Number of 3RDM elements
!> @param[in]  g3     3RDM
!> @param[in]  g2     2RDM
!> @param[in]  g1     1RDM
!> @param[in]  f3     contracted Fock matrix with 4RDM
!> @param[in]  f2     contracted Fock matrix with 3RDM
!> @param[in]  f1     contracted Fock matrix with 2RDM
!> @param[in]  iroot  MCSCF root number.
#ifdef _HDF5_
subroutine load_fciqmc_mats(idxG3,g3,g2,g1,f3,f2,f1,iroot,nLev)

  integer(kind=iwp), intent(in) :: iroot, nLev
  integer(kind=byte), intent(in) :: idxG3(6,nG3)
  real(kind=wp), intent(inout) :: g3(*), g2(nLev,nLev,nLev,nLev), g1(nLev,nLev), f3(*), f2(nLev,nLev,nLev,nLev), f1(nLev,nLev)
  integer(kind=iwp) :: i, t, u, v, x, y, z
  real(kind=wp), allocatable :: f3_temp(:,:,:,:,:,:), g3_temp(:,:,:,:,:,:)

  write(u6,'(a)') "Initiate transfer of CASPT2 intermediates."
  call mma_allocate(f3_temp,nLev,nLev,nLev,nLev,nLev,nLev,Label='f3_temp')
  call mma_allocate(g3_temp,nLev,nLev,nLev,nLev,nLev,nLev,Label='g3_temp')
  call load_six_tensor(g3_temp,'/spinfree/3300/',iroot,nLev)
  call load_six_tensor(f3_temp,'/spinfree/4400f/',iroot,nLev)

  if (TransformToNormalOrder) then
    ! \sum_{ab} f_{ab} e_{tu,vx,yz} E_{ab} -> \sum_{ab} f_{ab} e_{tu,vx,yz,ab}
    write(u6,'(a)') "Transform F.4RDM from factorised to normal order."
    call transform_f4rdm_normal_order(f3_temp,g3_temp)
  end if

  if (NonDiagonal) then
    write(u6,'(a)') "Transform intermediates to pseudo-canonical orbitals."
    call transform_six_index(g3_temp)
    call transform_six_index(f3_temp)
  end if

  write(u6,'(a)') "Trace lower rank intermediates from 3RDM and F.4RDM."
  call calc_f2_and_g2(f3_temp,g3_temp,f2,g2)
  call calc_f1_and_g1(f2,g2,f1,g1)

  ! convert into flattened arrays
  do i=1,nG3
    t = idxG3(1,i)
    u = idxG3(2,i)
    v = idxG3(3,i)
    x = idxG3(4,i)
    y = idxG3(5,i)
    z = idxG3(6,i)
    g3(i) = g3_temp(t,u,v,x,y,z)
    f3(i) = f3_temp(t,u,v,x,y,z)
  end do
  call mma_deallocate(f3_temp)
  call mma_deallocate(g3_temp)

contains

  subroutine transform_f4rdm_normal_order(f3,g3)
    ! The transformation should be
    ! \sum_{ab} f_{ab} e_{tu,vx,yz,ab} =
    !   \sum_{ab} f_{ab} [ e_{tu,vx,yz} E_{ab} \
    !                      - \delta_{az} e_{tu,vx,yb} \ (1)
    !                      - \delta_{ax} e_{tu,vb,yz} \ (2)
    !                      - \delta_{au} e_{tb,vx,yz} ].(3)

    real(kind=wp), intent(inout) :: f3(nLev,nLev,nLev,nLev,nLev,nLev)
    real(kind=wp), intent(in) :: g3(nLev,nLev,nLev,nLev,nLev,nLev)
    integer(kind=iwp) :: t, u, v, x, y, z
    real(kind=wp), allocatable :: fock_eigvecs(:,:), fockmat(:,:)

    call mma_allocate(fockmat,nLev,nLev,Label='fockmat')
    call mma_allocate(fock_eigvecs,nLev,nLev,Label='fock_eigvecs')
    call load_fockmat(fockmat,fock_eigvecs,nLev)
    do z=1,nLev
      do y=1,nLev
        do x=1,nLev
          do v=1,nLev
            do u=1,nLev
              do t=1,nLev
                f3(t,u,v,x,y,z) = f3(t,u,v,x,y,z)-sum(fockmat(z,:)*g3(t,u,v,x,y,:))-sum(fockmat(x,:)*g3(t,u,v,:,y,z))- &
                                  sum(fockmat(u,:)*g3(t,:,v,x,y,z))
                ! (1)
                !f3(t,u,v,x,y,z) = f3(t,u,v,x,y,z)-DDOT_(nLev,fockmat(z,:),1,g3(t,u,v,x,y,:),1)
                ! (2)
                !f3(t,u,v,x,y,z) = f3(t,u,v,x,y,z)-DDOT_(nLev,fockmat(x,:),1,g3(t,u,v,:,y,z),1)
                ! (3)
                !f3(t,u,v,x,y,z) = f3(t,u,v,x,y,z)-DDOT_(nLev,fockmat(u,:),1,g3(t,:,v,x,y,z),1)
              end do
            end do
          end do
        end do
      end do
    end do
    call mma_deallocate(fockmat)
    call mma_deallocate(fock_eigvecs)

  end subroutine

  !> @brief
  !>   Transform 3RDM and F.4RDM to pseudo-canonical orbitals.
  !>
  !> @param[inout]  g3    dense redundant 3RDM
  !> @param[inout]  f3    dense redundant F.4RDM
  !> @param[in]     nLev  number of levels
  subroutine transform_six_index(six_index)

    real(kind=wp), intent(inout) :: six_index(nLev,nLev,nLev,nLev,nLev,nLev)
    integer(kind=iwp) :: iter, u2, v2, x2, y2, z2  ! prevent sharing scope with upper function
    real(kind=wp), allocatable :: buffer(:), buffer2(:), fock_eigvecs(:,:), fockmat(:,:)

    call mma_allocate(fockmat,nLev,nLev,Label='fockmat')
    call mma_allocate(fock_eigvecs,nLev,nLev,Label='fock_eigvecs')
    call mma_allocate(buffer,nLev,Label='buffer')
    call mma_allocate(buffer2,nLev,Label='buffer2')
    call load_fockmat(fockmat,fock_eigvecs,nLev)
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
                    call dgemv_('T',nLev,nLev,One,fock_eigvecs,nLev,buffer,1,Zero,buffer2,1)
                    six_index(:,u2,v2,x2,y2,z2) = buffer2(:)
                  case (2)
                    buffer(:) = six_index(u2,:,v2,x2,y2,z2)
                    call dgemv_('T',nLev,nLev,One,fock_eigvecs,nLev,buffer,1,Zero,buffer2,1)
                    six_index(u2,:,v2,x2,y2,z2) = buffer2(:)
                  case (3)
                    buffer(:) = six_index(u2,v2,:,x2,y2,z2)
                    call dgemv_('T',nLev,nLev,One,fock_eigvecs,nLev,buffer,1,Zero,buffer2,1)
                    six_index(u2,v2,:,x2,y2,z2) = buffer2(:)
                  case (4)
                    buffer(:) = six_index(u2,v2,x2,:,y2,z2)
                    call dgemv_('T',nLev,nLev,One,fock_eigvecs,nLev,buffer,1,Zero,buffer2,1)
                    six_index(u2,v2,x2,:,y2,z2) = buffer2(:)
                  case (5)
                    buffer(:) = six_index(u2,v2,x2,y2,:,z2)
                    call dgemv_('T',nLev,nLev,One,fock_eigvecs,nLev,buffer,1,Zero,buffer2,1)
                    six_index(u2,v2,x2,y2,:,z2) = buffer2(:)
                  case (6)
                    buffer(:) = six_index(u2,v2,x2,y2,z2,:)
                    call dgemv_('T',nLev,nLev,One,fock_eigvecs,nLev,buffer,1,Zero,buffer2,1)
                    six_index(u2,v2,x2,y2,z2,:) = buffer2(:)
                end select
              end do
            end do
          end do
        end do
      end do
    end do
    call mma_deallocate(fockmat)
    call mma_deallocate(fock_eigvecs)
    call mma_deallocate(buffer)
    call mma_deallocate(buffer2)

  end subroutine transform_six_index

  pure subroutine calc_f2_and_g2(f3_temp,g3_temp,f2,g2)

    real(kind=wp), intent(in) :: f3_temp(nLev,nLev,nLev,nLev,nLev,nLev), g3_temp(nLev,nLev,nLev,nLev,nLev,nLev)
    real(kind=wp), intent(inout) :: f2(nLev,nLev,nLev,nLev), g2(nLev,nLev,nLev,nLev)
    integer(kind=iwp) :: t, u, v, x, w

    f2(:,:,:,:) = Zero
    g2(:,:,:,:) = Zero
    do w=1,nLev
      do x=1,nLev
        do v=1,nLev
          do u=1,nLev
            do t=1,nLev
              f2(t,u,v,x) = f2(t,u,v,x)+f3_temp(t,u,v,x,w,w)
              g2(t,u,v,x) = g2(t,u,v,x)+g3_temp(t,u,v,x,w,w)
            end do
          end do
        end do
      end do
    end do
    f2(:,:,:,:) = f2(:,:,:,:)/real(nActel-3,kind=wp)
    g2(:,:,:,:) = g2(:,:,:,:)/real(nActel-2,kind=wp)

  end subroutine calc_f2_and_g2

  pure subroutine calc_f1_and_g1(f2,g2,f1,g1)

    real(kind=wp), intent(in) :: f2(nLev,nLev,nLev,nLev), g2(nLev,nLev,nLev,nLev)
    real(kind=wp), intent(inout) :: f1(nLev,nLev), g1(nLev,nLev)
    integer(kind=iwp) :: t, u, w

    f1(:,:) = Zero
    g1(:,:) = Zero
    do w=1,nLev
      do u=1,nLev
        do t=1,nLev
          f1(t,u) = f1(t,u)+f2(t,u,w,w)
          g1(t,u) = g1(t,u)+g2(t,u,w,w)
        end do
      end do
    end do
    f1(:,:) = f1(:,:)/real(nActel-2,kind=wp)
    g1(:,:) = g1(:,:)/real(nActel-1,kind=wp)

  end subroutine calc_f1_and_g1

end subroutine load_fciqmc_mats

! required for MPI parallelisation
subroutine broadcast_filename(InFile)

  use filesystem, only: get_errno_, strerror_, symlink_

  character(len=*), intent(in) :: InFile
  character(len=1024) :: master
  integer(kind=iwp) :: lmaster1, err

  call prgmtranslate_master(InFile,master,lmaster1)
  call symlink_(trim(master),trim(InFile),err)
  if (err == 0) write(u6,*) strerror_(get_errno_())

end subroutine broadcast_filename

subroutine load_six_tensor(tensor,dataset,iroot,nLev)

  integer(kind=iwp), intent(in) :: iroot, nLev
  real(kind=wp), intent(inout) :: tensor(nLev,nLev,nLev,nLev,nLev,nLev)
  character(len=*), intent(in) :: dataset
  integer(kind=iwp) :: hdf5_dset, hdf5_file, hdf5_group, i, len6index(2), t, u, v, x, y, z
  logical(kind=iwp) :: tExist
  integer(kind=iwp), allocatable :: indices(:,:)
  real(kind=wp), allocatable :: values(:)

  if (myRank /= 0) call broadcast_filename('fciqmc.caspt2.'//str(iroot)//'.h5')
  call f_Inquire('fciqmc.caspt2.'//str(iroot)//'.h5',tExist)
  call verify_(tExist,'fciqmc.caspt2.'//str(iroot)//'.h5 does not exist.')
  hdf5_file = mh5_open_file_r('fciqmc.caspt2.'//str(iroot)//'.h5')
  hdf5_group = mh5_open_group(hdf5_file,trim(dataset))
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

  tensor(:,:,:,:,:,:) = Zero
  do i=1,len6index(2)
    t = indices(1,i)+1
    u = indices(4,i)+1
    v = indices(2,i)+1
    x = indices(5,i)+1
    y = indices(3,i)+1
    z = indices(6,i)+1
    tensor(t,u,v,x,y,z) = values(i)
    ! pre-contracted F4RDM is no longer hermitian
    if (dataset == '/spinfree/4400f/' .and. TransformToNormalOrder) then
      call apply_6fold_symmetry(tensor,t,u,v,x,y,z,values(i))
    else
      call apply_12fold_symmetry(tensor,t,u,v,x,y,z,values(i))
    end if
  end do

  call mma_deallocate(indices)
  call mma_deallocate(values)
  call mh5_close_file(hdf5_file)

contains

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

  pure subroutine apply_6fold_symmetry(array,t,u,v,x,y,z,val)
    ! If F4RDM is calculated from a histogrammed wave function,
    ! the 4RDM is rewritten as
    !     <e_pq,rs,tu E_vx>
    ! such that E_vx can be precontracted with |hist>. This change voids
    ! the hermiticity property.

    real(kind=wp), intent(inout) :: array(:,:,:,:,:,:)
    integer(kind=iwp), intent(in) :: t, u, v, x, y, z
    real(kind=wp), intent(in) :: val

    array(t,u,v,x,y,z) = val
    array(t,u,y,z,v,x) = val
    array(v,x,t,u,y,z) = val
    array(v,x,y,z,t,u) = val
    array(y,z,t,u,v,x) = val
    array(y,z,v,x,t,u) = val

  end subroutine apply_6fold_symmetry

end subroutine load_six_tensor

subroutine user_barrier()

  integer(kind=iwp) :: err
  logical(kind=iwp) :: proceed_found
  character(len=1024) :: WorkDir
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
    write(u6,'(4x,a)') 'Use the same FciDump as for the preceding CASCI.'
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
  if (myrank == 0) then
    write(u6,'(a)') 'PROCEED file found. Continuing with CASPT2.'
  else if (myRank /= 0) then
    call broadcast_filename('fciqmc.caspt2.'//str(mstate(jstate))//'.h5')
  end if

end subroutine user_barrier

subroutine load_fockmat(fock_matrix,fock_eigenvectors,nLev)
  ! sometimes eigenvectors are superfluous, but I/O should stay in one place
  ! and loading them is basically for free.

  integer(kind=iwp), intent(in) :: nLev
  real(kind=wp), intent(inout) :: fock_matrix(nLev,nLev), fock_eigenvectors(nLev,nLev)
  integer(kind=iwp) :: hdf5_dset, hdf5_file, hdf5_group, i, len2index(2), t, u
  logical(kind=iwp) :: tExist
  integer(kind=iwp), allocatable :: indices(:,:)
  real(kind=wp), allocatable :: values(:)

  if (myRank /= 0) call broadcast_filename('fockdump.h5')
  call f_Inquire('fockdump.h5',tExist)
  call verify_(tExist,'fockdump.h5 does not exist.')
  hdf5_file = mh5_open_file_r('fockdump.h5')
  hdf5_group = mh5_open_group(hdf5_file,'/')
  call mh5_fetch_dset(hdf5_group,'ACT_FOCK_EIGVECS',fock_eigenvectors)
  hdf5_dset = mh5_open_dset(hdf5_group,'ACT_FOCK_INDEX')
  len2index(:) = 0
  call mh5_get_dset_dims(hdf5_dset,len2index)
  call mma_allocate(indices,2,len2index(2))
  call mma_allocate(values,len2index(2))
  indices(:,:) = 0
  values(:) = Zero
  call mh5_fetch_dset(hdf5_group,'ACT_FOCK_VALUES',values)
  call mh5_fetch_dset(hdf5_group,'ACT_FOCK_INDEX',indices)
  call mh5_close_group(hdf5_group)
  call mh5_close_file(hdf5_file)

  fock_matrix(:,:) = Zero
  do i=1,len2index(2)
    t = indices(1,i)
    u = indices(2,i)  ! unlike M7 hdf5 files, fockdump is 1-base indexed
    fock_matrix(t,u) = values(i)
    fock_matrix(u,t) = values(i)
  end do
  call mma_deallocate(indices)
  call mma_deallocate(values)

end subroutine load_fockmat
#endif

end module fciqmc_interface
