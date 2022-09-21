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
! Copyright (C) 2020, Roland Lindh                                     *
!***********************************************************************

module Center_Info

use Definitions, only: iwp

implicit none
private

public :: Center_Info_Dmp, Center_Info_Free, Center_Info_Get, Center_Info_Init, dc, n_dc

#include "Molcas.fh"

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! D E C L A R E   D E R I V E D   T Y P E S
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

type Distinct_centers
  integer(kind=iwp) :: iChCnt = 0
  integer(kind=iwp) :: iStab(0:7) = 0
  integer(kind=iwp) :: nStab = 0
  integer(kind=iwp) :: iCoSet(0:7,0:7) = 0
  character(len=LenIn4) :: LblCnt = ''
end type Distinct_centers

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! E N D   D E C L A R E   D E R I V E D   T Y P E S
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

integer(kind=iwp), parameter :: nFields = 1+8+1+8**2
integer(kind=iwp) :: n_dc = 0
logical(kind=iwp) :: Initiated = .false.
type(Distinct_centers), allocatable :: dc(:)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Private extensions to mma interfaces

interface cptr2loff
  module procedure dc_cptr2loff
end interface
interface mma_Allocate
  module procedure dc_mma_allo_1D, dc_mma_allo_1D_lim
end interface
interface mma_Deallocate
  module procedure dc_mma_free_1D
end interface

contains

!***********************************************************************
!***********************************************************************
!
! This to make either the initial allocation of dc according to the default sizes
! as defined by the parameters in Molcas.fh or according to the actual sizes as recorded on the
! run file.

subroutine Center_Info_Init()

  use Definitions, only: u6

# include "macros.fh"
  unused_proc(mma_allocate(dc,[0,0]))

  if (Initiated) then
    write(u6,*) 'Center_Info already initiated!'
    write(u6,*) 'May the is a missing call to Center_Info_Free.'
    call Abend()
  end if
# ifdef _DEBUGPRINT_
  write(u6,*)
  write(u6,*) 'Enter Center_Info_Init'
# endif
  if (n_dc == 0) then
    call mma_allocate(dc,MxAtom,label='dc')
  else
    call mma_allocate(dc,n_dc,label='dc')
  end if
  Initiated = .true.
# ifdef _DEBUGPRINT_
  write(u6,*) 'Exit Center_Info_Init'
# endif

  return

end subroutine Center_Info_Init

!***********************************************************************
!***********************************************************************

subroutine Center_Info_Dmp()

  use stdalloc, only: mma_allocate, mma_deallocate
# ifdef _DEBUGPRINT_
  use Definitions, only: u6
# endif

  integer(kind=iwp) :: i, j, lcDmp, licDmp
  integer(kind=iwp), allocatable :: iDmp(:)
  character(len=LenIn4), allocatable :: cDmp(:)

  ! Integer dc stuff

# ifdef _DEBUGPRINT_
  write(u6,*)
  write(u6,*) 'Enter Center_Info_Dmp'
  write(u6,*) 'n_dc=',n_dc
# endif
  licDmp = n_dc*nFields
  call mma_allocate(iDmp,licDmp+1,Label='iDmp')
  j = 0
  do i=1,n_dc
    iDmp(j+1) = dc(i)%iChCnt
    j = j+1
    iDmp(j+1:j+8) = dc(i)%iStab(0:7)
    j = j+8
    iDmp(j+1) = dc(i)%nStab
    j = j+1
    iDmp(j+1:j+8*8) = reshape(dc(i)%iCoSet(:,:),[8*8])
    j = j+8*8
  end do
  iDmp(licDmp+1) = n_dc
  call Put_iArray('icDmp',iDmp,licDmp+1)
  call mma_deallocate(iDmp)

  call mma_allocate(cDmp,n_dc,Label='cDmp')
  do i=1,n_dc
    cDmp(i) = dc(i)%LblCnt
  end do
  lcDmp = n_dc*LenIn4
# ifdef _DEBUGPRINT_
  write(u6,*) 'cDmp=',cDmp(1:lcDmp)
# endif
  call Put_cArray('dc: cDmp',cDmp(1),lcDmp)
  call mma_deallocate(cDmp)
# ifdef _DEBUGPRINT_
  write(u6,*) 'Exit Center_Info_Dmp'
# endif

  return

end subroutine Center_Info_Dmp

!***********************************************************************
!***********************************************************************

subroutine Center_Info_Get()

  use stdalloc, only: mma_allocate, mma_deallocate
  use Definitions, only: u6

  integer(kind=iwp) :: i, j, lcDmp, Len1
  logical(kind=iwp) :: Found
  integer(kind=iwp), allocatable :: iDmp(:)
  character(len=LenIn4), allocatable :: cDmp(:)

# ifdef _DEBUGPRINT_
  write(u6,*)
  write(u6,*) 'Enter Center_Info_Get'
# endif
  call qpg_iArray('icDmp',Found,Len1)
  call mma_Allocate(iDmp,Len1,Label='iDmp')
  if (Found) then
    call Get_iArray('icDmp',iDmp,Len1)
  else
    write(u6,*) 'Center_Info_Get: icDmp not found!'
    call Abend()
  end if
  lcDmp = Len1-1
  n_dc = lcDmp/nFields

  ! Initiate the memory allocation of dc

  if (.not. Initiated) call Center_Info_Init()
# ifdef _DEBUGPRINT_
  write(u6,*) 'iDmp(1:Len1)=',iDmp(1:Len1)
  write(u6,*) 'Len1=',Len1
  write(u6,*) 'n_dc=',n_dc
  write(u6,*) 'lcDmp=',lcDmp
# endif
  j = 0
  do i=1,n_dc
    dc(i)%iChCnt = iDmp(j+1)
    j = j+1
    dc(i)%iStab(0:7) = iDmp(j+1:j+8)
    j = j+8
    dc(i)%nStab = iDmp(j+1)
    j = j+1
    dc(i)%iCoSet(:,:) = reshape(iDmp(j+1:j+8*8),[8,8])
    j = j+8*8
  end do
  call mma_deAllocate(iDmp)

  lcDmp = n_dc*LenIn4
# ifdef _DEBUGPRINT_
  write(u6,*) 'lcDmp=',lcDmp
# endif

  call qpg_cArray('dc: cDmp',Found,Len1)
  if (Len1 /= lcDmp) then
    write(u6,*) 'Center_Info_Get: Len1 /= lcDmp'
    call Abend()
  end if
  call mma_Allocate(cDmp,lcDmp,Label='cDmp')
  call Get_cArray('dc: cDmp',cDmp,lcDmp)
# ifdef _DEBUGPRINT_
  write(u6,*) 'cDmp=',cDmp(1:lcDmp)
# endif
  do i=1,n_dc
    dc(i)%LblCnt = cDmp(i)
  end do
  call mma_deAllocate(cDmp)
# ifdef _DEBUGPRINT_
  write(u6,*) 'Exit Center_Info_Get'
# endif

end subroutine Center_Info_Get

!***********************************************************************
!***********************************************************************

subroutine Center_Info_Free()

# ifdef _DEBUGPRINT_
  use Definitions, only: u6
# endif

  ! Deallocate all allocatable parts of dc.

  if (.not. allocated(dc)) return
# ifdef _DEBUGPRINT_
  write(u6,*)
  write(u6,*) 'Enter Center_Info_Free'
# endif
  call mma_deallocate(dc)
  n_dc = 0
  Initiated = .false.

# ifdef _DEBUGPRINT_
  write(u6,*) 'Exit Center_Info_Free'
# endif

  return

end subroutine Center_Info_Free

!***********************************************************************
!***********************************************************************

! Private extensions to mma_interfaces, using preprocessor templates
! (see src/mma_util/stdalloc.f)

! Define dc_cptr2loff, dc_mma_allo_1D, dc_mma_allo_1D_lim, dc_mma_free_1D
! (using _NO_GARBLE_ because all members are initialized)
#define _TYPE_ type(Distinct_centers)
#  define _NO_GARBLE_
#  define _FUNC_NAME_ dc_cptr2loff
#  include "cptr2loff_template.fh"
#  undef _FUNC_NAME_
#  define _SUBR_NAME_ dc_mma
#  define _DIMENSIONS_ 1
#  define _DEF_LABEL_ 'dc_mma'
#  include "mma_allo_template.fh"
#  undef _SUBR_NAME_
#  undef _DIMENSIONS_
#  undef _DEF_LABEL_
#  undef _NO_GARBLE_
#undef _TYPE_

end module Center_Info
