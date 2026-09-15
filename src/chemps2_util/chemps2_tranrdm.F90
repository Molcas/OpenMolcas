!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
!***********************************************************************
! Transform CheMPS2 RDMs from the optimized active-orbital basis to the
! pseudocanonical basis used internally by CASPT2.  The transformed data
! are stored in separate .tran files, preserving the original RDMs.

subroutine chemps2_tran2pdm(nact,xmat,chemroot)
#ifdef _MOLCAS_MPP_
use MPI, only: MPI_COMM_WORLD
use Para_Info, only: Is_Real_Par, King
use Definitions, only: MPIInt
#endif
use Definitions, only: wp, iwp, u6
implicit none
integer(kind=iwp), intent(in) :: nact, chemroot
real(kind=wp), intent(in) :: xmat(nact,nact)
#ifdef _MOLCAS_MPP_
integer(kind=MPIInt) :: mpi_error
if (King() .or. (.not. Is_Real_Par())) call transform_file(4,'2-RDM','molcas_2rdm.h5')
if (Is_Real_Par()) call MPI_Barrier(MPI_COMM_WORLD,mpi_error)
#else
call transform_file(4,'2-RDM','molcas_2rdm.h5')
#endif
contains
  subroutine transform_file(rank,group_name,base_name)
    use mh5, only: mh5_open_file_r, mh5_open_file_rw, mh5_open_group, mh5_open_dset, &
                  mh5_fetch_dset, mh5_put_dset, mh5_close_dset, mh5_close_group, mh5_close_file
    use stdalloc, only: mma_allocate, mma_deallocate
    integer(kind=iwp), intent(in) :: rank
    character(len=*), intent(in) :: group_name, base_name
    character(len=64) :: source_file, target_file
    character(len=10) :: rootindex
    integer(kind=iwp) :: file_h5, group_h5, dset_h5, ierr
    logical(kind=iwp) :: exists
    real(kind=wp), allocatable :: rdm(:)
    write(rootindex,'(i2)') chemroot-1
    source_file=trim(base_name)//'.r'//trim(adjustl(rootindex))
    target_file=trim(source_file)//'.tran'
    call f_inquire(trim(source_file),exists)
    if (.not. exists) then
      write(u6,*) 'CHEMPS2> Missing RDM file: ',trim(source_file)
      call abend()
    end if
    call mma_allocate(rdm,nact**rank,label='chemps2_rdm')
    file_h5=mh5_open_file_r(trim(source_file))
    group_h5=mh5_open_group(file_h5,group_name)
    call mh5_fetch_dset(group_h5,'elements',rdm)
    call mh5_close_group(group_h5)
    call mh5_close_file(file_h5)
    call transform_tensor(rank,rdm)
    call systemf('cp -f '//trim(source_file)//' '//trim(target_file),ierr)
    if (ierr /= 0) then
      write(u6,*) 'CHEMPS2> Could not create ',trim(target_file)
      call abend()
    end if
    file_h5=mh5_open_file_rw(trim(target_file))
    group_h5=mh5_open_group(file_h5,group_name)
    dset_h5=mh5_open_dset(group_h5,'elements')
    call mh5_put_dset(dset_h5,rdm)
    call mh5_close_dset(dset_h5)
    call mh5_close_group(group_h5)
    call mh5_close_file(file_h5)
    call mma_deallocate(rdm)
  end subroutine transform_file

  subroutine transform_tensor(rank,rdm)
    use stdalloc, only: mma_allocate, mma_deallocate
    integer(kind=iwp), intent(in) :: rank
    real(kind=wp), intent(inout) :: rdm(nact**rank)
    integer(kind=iwp) :: block, group, i, iter, offset
    real(kind=wp), allocatable :: out(:), vin(:), vout(:)
    block=nact**(rank-1)
    call mma_allocate(out,nact**rank,label='chemps2_rdm_out')
    call mma_allocate(vin,nact,label='chemps2_rdm_vin')
    call mma_allocate(vout,nact,label='chemps2_rdm_vout')
    do iter=1,rank
      do group=1,block
        offset=(group-1)*nact
        vin(:)=rdm(offset+1:offset+nact)
        call dgemv_('T',nact,nact,1.0_wp,xmat,nact,vin,1,0.0_wp,vout,1)
        do i=1,nact
          out(group+(i-1)*block)=vout(i)
        end do
      end do
      rdm(:)=out(:)
    end do
    call mma_deallocate(vout)
    call mma_deallocate(vin)
    call mma_deallocate(out)
  end subroutine transform_tensor
end subroutine chemps2_tran2pdm

subroutine chemps2_tran3pdm(nact,xmat,chemroot,tran3rdm)
#ifdef _MOLCAS_MPP_
use MPI, only: MPI_COMM_WORLD
use Para_Info, only: Is_Real_Par, King
use Definitions, only: MPIInt
#endif
use Definitions, only: wp, iwp, u6
implicit none
integer(kind=iwp), intent(in) :: nact, chemroot
real(kind=wp), intent(in) :: xmat(nact,nact)
logical(kind=iwp), intent(in) :: tran3rdm
character(len=8) :: group_name
character(len=18) :: base_name
#ifdef _MOLCAS_MPP_
integer(kind=MPIInt) :: mpi_error
#endif
if (tran3rdm) then
  group_name='3-RDM'
  base_name='molcas_3rdm.h5'
else
  group_name='F.4-RDM'
  base_name='molcas_f4rdm.h5'
end if
#ifdef _MOLCAS_MPP_
if (King() .or. (.not. Is_Real_Par())) call transform_file()
if (Is_Real_Par()) call MPI_Barrier(MPI_COMM_WORLD,mpi_error)
#else
call transform_file()
#endif
contains
  subroutine transform_file()
    use mh5, only: mh5_open_file_r, mh5_open_file_rw, mh5_open_group, mh5_open_dset, &
                  mh5_fetch_dset, mh5_put_dset, mh5_close_dset, mh5_close_group, mh5_close_file
    use stdalloc, only: mma_allocate, mma_deallocate
    character(len=64) :: source_file, target_file
    character(len=10) :: rootindex
    integer(kind=iwp) :: file_h5, group_h5, dset_h5, ierr
    logical(kind=iwp) :: exists
    real(kind=wp), allocatable :: rdm(:), out(:), vin(:), vout(:)
    integer(kind=iwp) :: block, group, i, iter, offset
    write(rootindex,'(i2)') chemroot-1
    source_file=trim(base_name)//'.r'//trim(adjustl(rootindex))
    target_file=trim(source_file)//'.tran'
    call f_inquire(trim(source_file),exists)
    if (.not. exists) then
      write(u6,*) 'CHEMPS2> Missing RDM file: ',trim(source_file)
      call abend()
    end if
    block=nact**5
    call mma_allocate(rdm,nact**6,label='chemps2_3rdm')
    call mma_allocate(out,nact**6,label='chemps2_3rdm_out')
    call mma_allocate(vin,nact,label='chemps2_3rdm_vin')
    call mma_allocate(vout,nact,label='chemps2_3rdm_vout')
    file_h5=mh5_open_file_r(trim(source_file))
    group_h5=mh5_open_group(file_h5,trim(group_name))
    call mh5_fetch_dset(group_h5,'elements',rdm)
    call mh5_close_group(group_h5)
    call mh5_close_file(file_h5)
    do iter=1,6
      do group=1,block
        offset=(group-1)*nact
        vin(:)=rdm(offset+1:offset+nact)
        call dgemv_('T',nact,nact,1.0_wp,xmat,nact,vin,1,0.0_wp,vout,1)
        do i=1,nact
          out(group+(i-1)*block)=vout(i)
        end do
      end do
      rdm(:)=out(:)
    end do
    call systemf('cp -f '//trim(source_file)//' '//trim(target_file),ierr)
    if (ierr /= 0) then
      write(u6,*) 'CHEMPS2> Could not create ',trim(target_file)
      call abend()
    end if
    file_h5=mh5_open_file_rw(trim(target_file))
    group_h5=mh5_open_group(file_h5,trim(group_name))
    dset_h5=mh5_open_dset(group_h5,'elements')
    call mh5_put_dset(dset_h5,rdm)
    call mh5_close_dset(dset_h5)
    call mh5_close_group(group_h5)
    call mh5_close_file(file_h5)
    call mma_deallocate(vout)
    call mma_deallocate(vin)
    call mma_deallocate(out)
    call mma_deallocate(rdm)
  end subroutine transform_file
end subroutine chemps2_tran3pdm
