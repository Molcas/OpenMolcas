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

use stdalloc, only: mma_allocate, mma_deallocate

implicit none

#ifdef _DMRG_
private

type qcm_names
  character(len=256), allocatable :: states(:) ! full checkpoint names for every state
end type

type(qcm_names), allocatable :: qcm_group_names(:)
character(len=256), allocatable :: qcm_prefixes(:)
! prefix for the particular group (although redundant but used in the new MPSSI interface)

public :: qcm_group_names
public :: qcm_prefixes
public :: qcmaquis_info_init, qcmaquis_info_deinit

! Private extension to mma interfaces

interface cptr2loff
  module procedure :: qcmn_cptr2loff
end interface
interface mma_allocate
  module procedure :: qcmn_mma_allo_1D, qcmn_mma_allo_1D_lim
end interface
interface mma_deallocate
  module procedure :: qcmn_mma_free_1D
end interface

contains

subroutine qcmaquis_info_init(igroup,nstates,tag)

  use Definitions, only: iwp, u6

  integer(kind=iwp), intent(in) :: igroup, nstates, tag
# ifdef _GARBLE_
  interface
    subroutine c_null_alloc(A)
      character(len=*), allocatable :: A(:)
    end subroutine c_null_alloc
  end interface
  integer(kind=iwp) :: i
# endif

  if (tag == 0) then
    call mma_allocate(qcm_group_names,igroup,label='qcm_group_names')
#   ifdef _GARBLE_
    ! Garbling corrupts the allocation status of allocatable components, use a hack to reset it
    do i=1,igroup
      call c_null_alloc(qcm_group_names(i)%states)
    end do
#   endif
    call mma_allocate(qcm_prefixes,igroup,label='qcm_prefixes')
  else if (tag == 1) then
    call mma_allocate(qcm_group_names(igroup)%states,nstates,label='qcm_igroup')
    qcm_group_names(igroup)%states = ''
  else if (tag == -1) then
    call mma_allocate(qcm_group_names,igroup,label='qcm_group_names')
#   ifdef _GARBLE_
    do i=1,igroup
      call c_null_alloc(qcm_group_names(i)%states)
    end do
#   endif
    call mma_allocate(qcm_group_names(igroup)%states,nstates,label='qcm_igroup')
    qcm_group_names(igroup)%states = ''
  else
    write(u6,*) 'unknown tag in qcmaquis_info_init'
    call abend()
  end if

# include "macros.fh"
  unused_proc(mma_allocate(qcm_group_names,[0,0]))

end subroutine qcmaquis_info_init

subroutine qcmaquis_info_deinit

  use Definitions, only: iwp

  integer(kind=iwp) :: i

  if (.not. allocated(qcm_group_names)) return
  do i=1,size(qcm_group_names)
    if (.not. allocated(qcm_group_names)) return
    if (allocated(qcm_group_names(i)%states)) call mma_deallocate(qcm_group_names(i)%states)
  end do
  call mma_deallocate(qcm_group_names)
  if(allocated(qcm_prefixes)) call mma_deallocate(qcm_prefixes)

end subroutine qcmaquis_info_deinit

! Private extension to mma_interfaces, using preprocessor templates
! (see src/mma_util/stdalloc.f)

! Define qcmn_cptr2loff, qcmn_mma_allo_1D, qcmn_mma_allo_1D_lim, qcmn_mma_free_1D
#define _TYPE_ type(qcm_names)
#  define _FUNC_NAME_ qcmn_cptr2loff
#  include "cptr2loff_template.fh"
#  undef _FUNC_NAME_
#  define _SUBR_NAME_ qcmn_mma
#  define _DIMENSIONS_ 1
#  define _DEF_LABEL_ 'qcm_mma'
#  include "mma_allo_template.fh"
#  undef _SUBR_NAME_
#  undef _DIMENSIONS_
#  undef _DEF_LABEL_
#undef _TYPE_
#endif
end module qcmaquis_info
