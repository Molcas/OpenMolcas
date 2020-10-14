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
! Copyright (C) 2017, Stefan Knecht                                    *
!***********************************************************************
module qcmaquis_info

 implicit none

 type qcm_names
 character(len=256), allocatable :: states(:) ! full checkpoint names for every state
 end type

 type(qcm_names), public, allocatable :: qcm_group_names(:)
 character(len=256), public, allocatable :: qcm_prefixes(:) ! prefix for the particular group (although redundant but used in the new MPSSI interface)

 save

 contains

 subroutine qcmaquis_info_init(igroup,nstates,tag)

   integer, intent(in) :: igroup, nstates, tag

   if(tag == 0)then
     allocate(qcm_group_names(igroup))
     allocate(qcm_prefixes(igroup))
     qcm_prefixes = ''

   else if(tag == 1)then
     allocate(qcm_group_names(igroup)%states(nstates)); qcm_group_names(igroup)%states = ''
   else if(tag == -1)then
     allocate(qcm_group_names(igroup))
     allocate(qcm_group_names(igroup)%states(nstates)); qcm_group_names(igroup)%states = ''
   else
     write(6,*) 'unknown tag in qcmaquis_info_init'
     call abend()
   end if
 end subroutine qcmaquis_info_init

 subroutine qcmaquis_info_deinit
   integer             :: i
   if(.not.allocated(qcm_group_names)) return
   do i = 1, size(qcm_group_names)
     if(.not.allocated(qcm_group_names)) return
     if(allocated(qcm_group_names(i)%states)) deallocate(qcm_group_names(i)%states)
   end do
   deallocate(qcm_group_names)
   if(allocated(qcm_prefixes)) deallocate(qcm_prefixes)
 end subroutine qcmaquis_info_deinit

end module qcmaquis_info
