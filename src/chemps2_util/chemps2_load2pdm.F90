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

subroutine chemps2_load2pdm(NAC,PT,CHEMROOT)

use mh5, only: mh5_open_file_r, mh5_open_group, mh5_fetch_dset, mh5_close_group, mh5_close_file
use stdalloc, only: mma_allocate, mma_deallocate
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: NAC, CHEMROOT
real(kind=wp), intent(out) :: PT(NAC,NAC,NAC,NAC)
character(len=30) :: file_2rdm
character(len=10) :: rootindex
integer(kind=iwp) :: file_h5, group_h5, i, j, k, l, idx
logical(kind=iwp) :: irdm
real(kind=wp), allocatable :: two_rdm(:)

write(rootindex,'(i2)') chemroot-1
file_2rdm = 'molcas_2rdm.h5.r'//trim(adjustl(rootindex))
file_2rdm = trim(adjustl(file_2rdm))
call f_inquire(file_2rdm,irdm)
if (.not. irdm) then
  write(u6,'(1x,a15,i3,a16)') 'CHEMPS2> Root: ',CHEMROOT,' :: No 2RDM file'
  call abend()
end if

!#ifdef _MOLCAS_MPP_
!if (MPP() .and. KING()) then
!#endif

call mma_allocate(two_rdm,NAC**4,label='two_rdm')

file_h5 = mh5_open_file_r(file_2rdm)
group_h5 = mh5_open_group(file_h5,'2-RDM')
call mh5_fetch_dset(group_h5,'elements',two_rdm)
call mh5_close_group(group_h5)
call mh5_close_file(file_h5)
!#ifdef _MOLCAS_MPP_
!end if
!call MPI_Bcast(two_rdm,NAC**4,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERROR4)
!#endif

do i=1,NAC
  do j=1,NAC
    do k=1,NAC
      do l=1,NAC
        idx = i+NAC*(k-1+NAC*(j-1+NAC*(l-1)))
        PT(i,j,k,l) = two_rdm(idx)
      end do
    end do
  end do
end do

call mma_deallocate(two_rdm)

return

end subroutine chemps2_load2pdm
