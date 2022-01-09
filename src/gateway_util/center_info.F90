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

implicit none
private

public :: dc, n_dc, Center_Info_Init, Center_Info_Dmp, Center_Info_Get, Center_Info_Free

#include "Molcas.fh"
#include "stdalloc.fh"

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! D E C L A R E   D E R I V E D   T Y P E S
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

type Distinct_centers
  !sequence
  integer :: iChCnt = 0
  integer :: iStab(0:7) = [0,0,0,0,0,0,0,0]
  integer :: nStab = 0
  integer :: iCoSet(0:7,0:7) = reshape([0,0,0,0,0,0,0,0, &
                                        0,0,0,0,0,0,0,0, &
                                        0,0,0,0,0,0,0,0, &
                                        0,0,0,0,0,0,0,0, &
                                        0,0,0,0,0,0,0,0, &
                                        0,0,0,0,0,0,0,0, &
                                        0,0,0,0,0,0,0,0, &
                                        0,0,0,0,0,0,0,0],[8,8])
  character(LEN=LENIN4) :: LblCnt = ''
end type Distinct_centers

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! E N D   D E C L A R E   D E R I V E D   T Y P E S
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

integer, parameter :: nFields = 1+8+1+64
logical :: Initiated = .false.
integer :: n_dc = 0
type(Distinct_centers), allocatable :: dc(:)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

interface
  subroutine Abend()
  end subroutine Abend
  subroutine Put_iArray(Label,data,nData)
    character*(*) Label
    integer nData
    integer data(nData)
  end subroutine Put_iArray
  subroutine Get_iArray(Label,data,nData)
    character*(*) Label
    integer nData
    integer data(nData)
  end subroutine Get_iArray
  subroutine Qpg_iArray(Label,Found,nData)
    character*(*) Label
    logical Found
    integer nData
  end subroutine Qpg_iArray
end interface

!***********************************************************************
!***********************************************************************

contains

!***********************************************************************
!***********************************************************************
!
! This to make either the initial allocation of dbsc and Shells according to the default sizes
! as defined by the parameters in Molcas.fh or according to the actual sizes as recorded on the
! run file.

subroutine Center_Info_Init()

  if (Initiated) then
    write(6,*) 'Center_Info already initiated!'
    write(6,*) 'May the is a missing call to Center_Info_Free.'
    call Abend()
  end if
# ifdef _DEBUGPRINT_
  write(6,*)
  write(6,*) 'Enter Center_Info_Init'
# endif
  if (n_dc == 0) then
    allocate(dc(1:MxAtom))
  else
    allocate(dc(1:n_dc))
  end if
  Initiated = .true.
# ifdef _DEBUGPRINT_
  write(6,*) 'Exit Center_Info_Init'
# endif

  return

end subroutine Center_Info_Init

!***********************************************************************
!***********************************************************************

subroutine Center_Info_Dmp()

  integer i, j, k, licDmp, lcDmp
  integer, allocatable :: iDmp(:)
  character(LEN=1), allocatable :: cDmp(:)

  ! Integer dc stuff

# ifdef _DEBUGPRINT_
  write(6,*)
  write(6,*) 'Enter Center_Info_Dmp'
  write(6,*) 'n_dc=',n_dc
# endif
  licDmp = n_dc*nFields
  call mma_Allocate(iDmp,licDmp+1,Label='iDmp')
  j = 0
  do i=1,n_dc
    iDmp(j+1) = dc(i)%iChCnt
    j = j+1
    iDmp(j+1:j+8) = dc(i)%iStab(0:7)
    j = j+8
    iDmp(j+1) = dc(i)%nStab
    j = j+1
    iDmp(j+1:j+8) = dc(i)%iCoSet(:,0)
    j = j+8
    iDmp(j+1:j+8) = dc(i)%iCoSet(:,1)
    j = j+8
    iDmp(j+1:j+8) = dc(i)%iCoSet(:,2)
    j = j+8
    iDmp(j+1:j+8) = dc(i)%iCoSet(:,3)
    j = j+8
    iDmp(j+1:j+8) = dc(i)%iCoSet(:,4)
    j = j+8
    iDmp(j+1:j+8) = dc(i)%iCoSet(:,5)
    j = j+8
    iDmp(j+1:j+8) = dc(i)%iCoSet(:,6)
    j = j+8
    iDmp(j+1:j+8) = dc(i)%iCoSet(:,7)
    j = j+8
  end do
  iDmp(licDmp+1) = n_dc
  call Put_iArray('icDmp',iDmp,licDmp+1)
  call mma_deallocate(iDmp)

  lcDmp = n_dc*LENIN4
# ifdef _DEBUGPRINT_
  write(6,*) 'lcDmp=',n_dc
# endif
  call mma_allocate(cDmp,lcDmp,Label='cDmp')
  k = 0
  do i=1,n_dc
    do j=1,LENIN4
      cDmp(k+j) = dc(i)%LblCnt(j:j)
    end do
    k = k+LENIN4
  end do
# ifdef _DEBUGPRINT_
  write(6,*) 'cDmp=',cDmp(1:lcDmp)
# endif
  call Put_cArray('dc: cDmp',cDmp(1),lcDmp)
  call mma_deallocate(cDmp)
# ifdef _DEBUGPRINT_
  write(6,*) 'Exit Center_Info_Dmp'
# endif

  return

end subroutine Center_Info_Dmp

!***********************************************************************
!***********************************************************************

subroutine Center_Info_Get()

  integer, allocatable :: iDmp(:)
  character(LEN=1), allocatable :: cDmp(:)
  logical Found
  integer i, j, k, Len, lcDmp

# ifdef _DEBUGPRINT_
  write(6,*)
  write(6,*) 'Enter Center_Info_Get'
# endif
  call qpg_iArray('icDmp',Found,Len)
  call mma_Allocate(iDmp,Len,Label='iDmp')
  if (Found) then
    call Get_iArray('icDmp',iDmp,Len)
  else
    write(6,*) 'Center_Info_Get: icDmp not found!'
    call Abend()
  end if
  lcDmp = Len-1
  n_dc = lcDmp/nFields

  ! Initiate the memory allocation of dc

  if (.not. Initiated) call Center_Info_Init()
# ifdef _DEBUGPRINT_
  write(6,*) 'iDmp(1:Len)=',iDmp(1:Len)
  write(6,*) 'Len=',Len
  write(6,*) 'n_dc=',n_dc
  write(6,*) 'lcDmp=',lcDmp
# endif
  j = 0
  do i=1,n_dc
    dc(i)%iChCnt = iDmp(j+1)
    j = j+1
    dc(i)%iStab(0:7) = iDmp(j+1:j+8)
    j = j+8
    dc(i)%nStab = iDmp(j+1)
    j = j+1
    dc(i)%iCoSet(:,0) = iDmp(j+1:j+8)
    j = j+8
    dc(i)%iCoSet(:,1) = iDmp(j+1:j+8)
    j = j+8
    dc(i)%iCoSet(:,2) = iDmp(j+1:j+8)
    j = j+8
    dc(i)%iCoSet(:,3) = iDmp(j+1:j+8)
    j = j+8
    dc(i)%iCoSet(:,4) = iDmp(j+1:j+8)
    j = j+8
    dc(i)%iCoSet(:,5) = iDmp(j+1:j+8)
    j = j+8
    dc(i)%iCoSet(:,6) = iDmp(j+1:j+8)
    j = j+8
    dc(i)%iCoSet(:,7) = iDmp(j+1:j+8)
    j = j+8
  end do
  call mma_deAllocate(iDmp)

  lcDmp = n_dc*LENIN4
# ifdef _DEBUGPRINT_
  write(6,*) 'lcDmp=',lcDmp
# endif

  call qpg_cArray('dc: cDmp',Found,Len)
  if (Len /= lcDmp) then
    write(6,*) 'Center_Info_Get: Len /= lcDmp'
    call Abend()
  end if
  call mma_Allocate(cDmp,lcDmp,Label='cDmp')
  call Get_cArray('dc: cDmp',cDmp(1),lcDmp)
# ifdef _DEBUGPRINT_
  write(6,*) 'cDmp=',cDmp(1:lcDmp)
# endif
  k = 0
  do i=1,n_dc
    do j=1,LENIN4
      dc(i)%LblCnt(j:j) = cDmp(k+j)
    end do
    k = k+LENIN4
  end do
  call mma_deAllocate(cDmp)
# ifdef _DEBUGPRINT_
  write(6,*) 'Exit Center_Info_Get'
# endif

end subroutine Center_Info_Get

!***********************************************************************
!***********************************************************************

subroutine Center_Info_Free()

  ! Deallocate all allocatable parts of dc.

  if (.not. allocated(dc)) return
# ifdef _DEBUGPRINT_
  write(6,*)
  write(6,*) 'Enter Center_Info_Free'
# endif
  deallocate(dc)
  n_dc = 0
  Initiated = .false.

# ifdef _DEBUGPRINT_
  write(6,*) 'Exit Center_Info_Free'
# endif

  return

end subroutine Center_Info_Free

!***********************************************************************
!***********************************************************************

end module Center_Info
