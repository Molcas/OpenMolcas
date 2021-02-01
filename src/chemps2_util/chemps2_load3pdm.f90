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
! Subroutine to load 3-RDM and F4-RDM
! Written by Quan Phung and Sebastian Wouters, Leuven, Aug 2016

subroutine chemps2_load3pdm( NAC, idxG3, NG3, storage, doG3, EPSA, F2, chemroot )

  USE ISO_C_BINDING
#ifdef _MOLCAS_MPP_
  USE MPI
#endif
  USE mh5, ONLY: mh5_open_file_r, mh5_open_group, mh5_fetch_dset, mh5_close_group, mh5_close_file

  IMPLICIT NONE
  INTEGER, INTENT(IN)   :: NAC, NG3, chemroot
  INTEGER*1, INTENT(IN) :: idxG3( 6, NG3 )
  REAL*8, INTENT(OUT)   :: storage(*)
  LOGICAL, INTENT(IN)   :: doG3
  REAL*8, INTENT(IN)    :: EPSA( NAC )
  REAL*8, INTENT(OUT)   :: F2 ( NAC, NAC, NAC, NAC )

  CHARACTER(LEN=30) :: file_3rdm
  CHARACTER(LEN=30) :: file_f4rdm
  LOGICAL           :: irdm, jrdm

  INTEGER            :: file_h5, group_h5
  character(len=10) :: rootindex

  INTEGER :: ip1, ip2, ip3, iq1, iq2, iq3, idx, iG3

  REAL*8, DIMENSION( NAC * NAC * NAC * NAC * NAC * NAC ) :: buffer

  write(rootindex,"(I2)") chemroot-1
  file_3rdm="molcas_3rdm.h5.r"//trim(adjustl(rootindex))
  file_f4rdm="molcas_f4rdm.h5.r"//trim(adjustl(rootindex))
  file_3rdm=trim(adjustl(file_3rdm))
  file_f4rdm=trim(adjustl(file_f4rdm))
  call f_inquire(file_3rdm, irdm)
  call f_inquire(file_f4rdm, jrdm)
  if ((.NOT. irdm) .OR. (.NOT. jrdm)) then
     write(6,'(1X,A15,I3,A26)') 'CHEMPS2> Root: ',CHEMROOT,' :: No 3RDM or F.4RDM file'
     call abend()
  endif

!#ifdef _MOLCAS_MPP_
!if ( MPP().AND.KING() ) then
!#endif
  If (doG3.EQV..true.) Then
    file_h5 = mh5_open_file_r(file_3rdm)
    group_h5 = mh5_open_group(file_h5, '3-RDM')
  Else
    file_h5 = mh5_open_file_r(file_f4rdm)
    group_h5 = mh5_open_group(file_h5, 'F.4-RDM')
  End If
  call mh5_fetch_dset(group_h5, 'elements', buffer)
  call mh5_close_group(group_h5)
  call mh5_close_file(file_h5)
!#ifdef _MOLCAS_MPP_
!end if
!call MPI_Bcast( buffer, NAC * NAC * NAC * NAC * NAC * NAC, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, IERROR4 )
!#endif

  do iG3=1,NG3
    ip1 = idxG3( 1, iG3 ) - 1
    iq1 = idxG3( 2, iG3 ) - 1
    ip2 = idxG3( 3, iG3 ) - 1
    iq2 = idxG3( 4, iG3 ) - 1
    ip3 = idxG3( 5, iG3 ) - 1
    iq3 = idxG3( 6, iG3 ) - 1
    idx = ip1 + NAC * ( ip2 + NAC * ( ip3 + NAC * ( iq1 + NAC * ( iq2 + NAC * iq3 ) ) ) )
    storage( iG3 ) = buffer( 1 + idx )
  end do

  if (doG3.EQV..true.) Then
    do iq2=1,NAC
      do ip2=1,NAC
        do iq1=1,NAC
          do ip1=1,NAC
            F2( ip1, iq1, ip2, iq2 ) = 0.0
            do ip3=1,NAC
              idx = ip1 + NAC * ( ip2 - 1 + NAC * ( ip3 - 1 + NAC * ( iq1 - 1 + NAC * ( iq2 - 1 + NAC * ( ip3 - 1 ) ) ) )  )
              F2( ip1, iq1, ip2, iq2 ) = F2( ip1, iq1, ip2, iq2 ) + EPSA( ip3 ) * buffer( idx )
            end do
          end do
        end do
      end do
    end do
  end if
end subroutine
