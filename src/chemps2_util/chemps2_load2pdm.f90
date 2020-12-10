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

  USE ISO_C_BINDING
#ifdef _MOLCAS_MPP_
  USE MPI
#endif

  IMPLICIT NONE
  INTEGER, INTENT(IN) :: NAC, CHEMROOT
  REAL*8, INTENT(OUT) :: PT( NAC, NAC, NAC, NAC )

  CHARACTER(LEN=30) :: file_2rdm

#include "mh5.fh"
  INTEGER            :: file_h5, group_h5
  LOGICAL            :: irdm

  INTEGER :: i,j,k,l,idx
  character(len=10) :: rootindex
#ifdef _MOLCAS_MPP_
  EXTERNAL Is_Real_Par, KING
  Logical KING
  Logical Is_Real_Par
#endif

  REAL*8, DIMENSION( NAC * NAC * NAC * NAC ) :: two_rdm

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

  file_h5 = mh5_open_file_r(file_2rdm)
  group_h5 = mh5_open_group(file_h5, '2-RDM')
  call mh5_fetch_dset_array_real(group_h5, 'elements', two_rdm)
  call mh5_close_group(group_h5)
  call mh5_close_file(file_h5)
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
