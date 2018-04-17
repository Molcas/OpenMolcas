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
! Copyright (C) 2016, Sebastian Wouters                                *
!               2016, Quan Phung                                       *
!***********************************************************************
! Subroutine to load 2RDM
! Written by Quan Phung and Sebastian Wouters, Leuven, Aug 2016

subroutine chemps2_load2pdm( NAC, PT, CHEMROOT )

  USE HDF5
  USE ISO_C_BINDING
#ifdef _MOLCAS_MPP_
  USE MPI
#endif

  IMPLICIT NONE
  INTEGER, INTENT(IN) :: NAC, CHEMROOT
  REAL*8, INTENT(OUT) :: PT( NAC, NAC, NAC, NAC )

  CHARACTER(LEN=30) :: file_2rdm

  INTEGER( HID_T )   :: file_h5, group_h5, space_h5, dset_h5 ! Handles
  INTEGER(4)         :: hdferr
  TYPE( C_PTR )      :: f_ptr
  LOGICAL            :: irdm

  INTEGER :: i,j,k,l,idx
  character(len=10) :: rootindex
#ifdef _MOLCAS_MPP_
  EXTERNAL Is_Real_Par, KING
  Logical KING
  Logical Is_Real_Par
#endif

  REAL*8, DIMENSION( 1 : NAC * NAC * NAC * NAC ), TARGET :: two_rdm

  write(rootindex,"(I2)") chemroot-1
  file_2rdm="molcas_2rdm.h5.r"//trim(adjustl(rootindex))
  file_2rdm=trim(adjustl(file_2rdm))
  call f_inquire(file_2rdm, irdm)
  if (.NOT. irdm) then
     write(6,'(1X,A15,I3,A16)') 'CHEMPS2> Root: ',CHEMROOT,' :: No 2RDM file'
     call abend()
  endif


!#ifdef _MOLCAS_MPP_
!  if ( MPP().AND.KING() ) then
!#endif

  CALL h5open_f( hdferr )
  CALL h5fopen_f( file_2rdm, H5F_ACC_RDONLY_F, file_h5, hdferr )
  CALL h5gopen_f( file_h5, "2-RDM", group_h5, hdferr )
  CALL h5dopen_f( group_h5, "elements", dset_h5, hdferr )
  CALL h5dget_space_f( dset_h5, space_h5, hdferr )
  f_ptr = C_LOC( two_rdm( 1 ) )
  CALL h5dread_f( dset_h5, H5T_NATIVE_DOUBLE, f_ptr, hdferr )
  CALL h5dclose_f( dset_h5 , hdferr )
  CALL h5sclose_f( space_h5, hdferr )
  CALL h5gclose_f( group_h5, hdferr )
  CALL h5fclose_f( file_h5,  hdferr )
!#ifdef _MOLCAS_MPP_
!  end if
!  call MPI_Bcast( two_rdm, NAC * NAC * NAC * NAC, MPI_DOUBLE_PRECISION, 0,     MPI_COMM_WORLD, IERROR4 )
!#endif


  do i=1,NAC
     do j=1,NAC
       do k=1,NAC
         do l=1,NAC
           idx = i + NAC * ( k - 1 + NAC * ( j - 1 + NAC * ( l - 1 )))
           PT( i, j, k, l ) = two_rdm( idx )
         enddo
       enddo
     enddo
  enddo

end subroutine
