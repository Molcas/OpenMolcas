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
! Copyright (C) 2018, Quan Phung                                       *
!***********************************************************************
! Load all elements of a CheMPS2 3-RDM. The flattened ordering is the
! native CheMPS2 ordering: p1,p2,p3,q1,q2,q3.

subroutine chemps2_load3pdm_all(NAC,storage,chemroot,TRANS)

use mh5, only: mh5_open_file_r, mh5_open_group, mh5_fetch_dset, mh5_close_group, mh5_close_file
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: NAC, chemroot
real(kind=wp), intent(out) :: storage(NAC**6)
logical(kind=iwp), intent(in) :: TRANS
character(len=64) :: file_3rdm
character(len=10) :: rootindex
integer(kind=iwp) :: file_h5, group_h5
logical(kind=iwp) :: rdm_found

write(rootindex,'(i2)') chemroot-1
file_3rdm = 'molcas_3rdm.h5.r'//trim(adjustl(rootindex))
if (TRANS) file_3rdm = trim(file_3rdm)//'.tran'

call f_inquire(trim(file_3rdm),rdm_found)
if (.not. rdm_found) then
  write(u6,'(1x,a15,i3,a17)') 'CHEMPS2> Root: ',chemroot,' :: No 3-RDM file'
  call abend()
end if

file_h5 = mh5_open_file_r(trim(file_3rdm))
group_h5 = mh5_open_group(file_h5,'3-RDM')
call mh5_fetch_dset(group_h5,'elements',storage)
call mh5_close_group(group_h5)
call mh5_close_file(file_h5)

end subroutine chemps2_load3pdm_all
