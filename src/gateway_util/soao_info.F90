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

module SOAO_Info

use Definitions, only: iwp

implicit none
private

integer(kind=iwp), allocatable :: iAOtSO(:,:), iSOInf(:,:)
integer(kind=iwp) :: iOffSO(0:7) = 0, nIrrep = 0, nSOInf = 0

! Do not change this value, since it also explicitly signals symmetry information.
! This is explicitly used in some routines. (where?)
integer(kind=iwp) :: initValue = -99999999

public :: iAOtSO, iOffSO, iSOInf, nSOInf, SOAO_Info_Dmp, SOAO_Info_Free, SOAO_Info_Get, SOAO_Info_Init

!***********************************************************************
!***********************************************************************

contains

!***********************************************************************
!***********************************************************************

subroutine SOAO_Info_Init(nSize,mIrrep)

  use stdalloc, only: mma_allocate

  integer(kind=iwp), intent(in) :: mIrrep, nSize

  if (allocated(iSOInf) .or. allocated(iAOtSO)) call SOAO_Info_Free()
  nSOInf = nSize
  nIrrep = mIrrep
  call mma_allocate(iSOInf,3,nSOInf,Label='iSOInf')
  iSOInf(:,:) = initValue
  call mma_allocate(iAOtSO,[1,nSOInf],[0,nIrrep-1],Label='iAOtSO')
  iAOtSO(:,:) = initValue

end subroutine SOAO_Info_Init

!***********************************************************************

subroutine SOAO_Info_Dmp()

  use stdalloc, only: mma_allocate, mma_deallocate
# ifdef _DEBUGPRINT_
  use Definitions, only: u6
# endif

  integer(kind=iwp) :: i, j
  integer(kind=iwp), allocatable :: iDmp(:)

# ifdef _DEBUGPRINT_
  write(u6,*)
  write(u6,*) 'Enter SOAO_Info_Dmp'
  write(u6,*)
  write(u6,*)
  write(u6,*) 'iAOtSO:'
  do i=0,nIrrep-1
    write(u6,'(8I9)') iAOtSO(1:nSOInf,i)
  end do
# endif
  call mma_allocate(iDmp,3*nSOInf+8,Label='iDmp')
  j = 0
  do i=1,nSOInf
    iDmp(j+1:j+3) = iSOInf(1:3,i)
    j = j+3
  end do
  iDmp(j+1:j+8) = iOffSO(0:7)
  call Put_iArray('iSOInf',iDmp,3*nSOInf+8)
  call mma_deallocate(iDmp)
  call Put_iArray('iAOtSO',iAOtSO,nSOInf*nIrrep)

end subroutine SOAO_Info_Dmp

!***********************************************************************

subroutine SOAO_Info_Get()

  use Definitions, only: u6

  use stdalloc, only: mma_allocate, mma_deallocate

  integer(kind=iwp) :: i, j
  logical(kind=iwp) :: Found
  integer(kind=iwp), allocatable :: iDmp(:)

# ifdef _DEBUGPRINT_
  write(u6,*)
  write(u6,*) 'Enter SOAO_Info_Get'
  write(u6,*)
# endif
  if (allocated(iSOInf) .or. allocated(iAOtSO)) call SOAO_Info_Free()
  call Qpg_iArray('iSOInf',Found,nSOInf)
  if (.not. Found) then
    write(u6,*) 'SOAO_Info_Get: iSOInf not found.'
    call Abend()
  end if
  nSOInf = (nSOInf-8)/3
  call mma_allocate(iSOInf,3,nSOInf,Label='iSOInf')
  call mma_allocate(iDmp,3*nSOInf+8,Label='iDmp')
  call Get_iArray('iSOInf',iDmp,3*nSOInf+8)
  j = 0
  do i=1,nSOInf
    iSOInf(1:3,i) = iDmp(j+1:j+3)
    j = j+3
  end do
  iOffSO(0:7) = iDmp(j+1:j+8)
  call mma_deallocate(iDmp)

  call Qpg_iArray('iAOtSO',Found,nIrrep)
  if (.not. Found) then
    write(u6,*) 'SOAO_Info_Get: iAOtSO not found.'
    call Abend()
  end if
  nIrrep = nIrrep/nSOInf
  call mma_allocate(iAOtSO,[1,nSOInf],[0,nIrrep-1],Label='iAOtSO')
  call Get_iArray('iAOtSO',iAOtSO,nSOInf*nIrrep)
# ifdef _DEBUGPRINT_
  write(u6,*)
  write(u6,*) 'iAOtSO:'
  do i=0,nIrrep-1
    write(u6,'(8I9)') iAOtSO(1:nSOInf,i)
  end do
  write(u6,*)
  write(u6,*) 'Exit SOAO_Info_Get'
  write(u6,*)
# endif

end subroutine SOAO_Info_Get

!***********************************************************************

subroutine SOAO_Info_Free()

  use stdalloc, only: mma_deallocate

  if (allocated(iSOInf)) call mma_deallocate(iSOInf)
  if (allocated(iAOtSO)) call mma_deallocate(iAOtSO)
  nSOInf = 0
  nIrrep = 0

end subroutine SOAO_Info_Free

!***********************************************************************
!***********************************************************************

end module SOAO_Info
