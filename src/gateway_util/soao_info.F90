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

implicit none
private

public :: iSOInf, iAOtSO, nSOInf, SOAO_Info_Init, SOAO_Info_Dmp, SOAO_Info_Get, SOAO_Info_Free, iOffSO

#include "stdalloc.fh"
#include "itmax.fh"
integer, allocatable :: iSOInf(:,:)
integer, allocatable :: iAOtSO(:,:)
integer :: iOffSO(0:7) = [0,0,0,0,0,0,0,0]
integer :: nSOInf = 0
integer :: nIrrep = 0

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

! This to make either the initial allocation of dbsc and Shells according to the default sizes
! as defined by the parameters in Molcas.fh or according to the actual sizes as recorded on the
! run file.

!***********************************************************************

subroutine SOAO_Info_Init(nSize,mIrrep)

  implicit none
  integer nSize, mIrrep

  if (allocated(iSOInf) .or. allocated(iAOtSO)) call SOAO_Info_Free()
  nSOInf = nSize
  nIrrep = mIrrep
  call mma_allocate(iSOInf,3,nSOInf,Label='iSOInf')
  iSOInf(:,:) = -99999999
  ! Do not change this value, since it also explicitly signals symmetry information.
  ! This is explicitly used in some routines.
  call mma_allocate(iAOtSO,[1,nSOInf],[0,nIrrep-1],Label='iAOtSO')
  iAOtSO(:,:) = -99999999        ! Dito

end subroutine SOAO_Info_Init

!***********************************************************************

subroutine SOAO_Info_Dmp()

  integer, allocatable :: iDmp(:)
  integer i, j

# ifdef _DEBUGPRINT_
  write(6,*)
  write(6,*) 'Enter SOAO_Info_Dmp'
  write(6,*)
  write(6,*)
  write(6,*) 'iAOtSO:'
  do i=0,nIrrep-1
    write(6,'(8I9)') iAOtSO(1:nSOInf,i)
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

  integer, allocatable :: iDmp(:)
  integer i, j
  logical Found

# ifdef _DEBUGPRINT_
  write(6,*)
  write(6,*) 'Enter SOAO_Info_Get'
  write(6,*)
# endif
  if (allocated(iSOInf) .or. allocated(iAOtSO)) call SOAO_Info_Free()
  call Qpg_iArray('iSOInf',Found,nSOInf)
  if (.not. Found) then
    write(6,*) 'SOAO_Info_Get: iSOInf not found.'
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
    write(6,*) 'SOAO_Info_Get: iAOtSO not found.'
    call Abend()
  end if
  nIrrep = nIrrep/nSOInf
  call mma_allocate(iAOtSO,[1,nSOInf],[0,nIrrep-1],Label='iAOtSO')
  call Get_iArray('iAOtSO',iAOtSO,nSOInf*nIrrep)
# ifdef _DEBUGPRINT_
  write(6,*)
  write(6,*) 'iAOtSO:'
  do i=0,nIrrep-1
    write(6,'(8I9)') iAOtSO(1:nSOInf,i)
  end do
  write(6,*)
  write(6,*) 'Exit SOAO_Info_Get'
  write(6,*)
# endif

end subroutine SOAO_Info_Get

!***********************************************************************

subroutine SOAO_Info_Free()

  if (allocated(iSOInf)) call mma_deallocate(iSOInf)
  if (allocated(iAOtSO)) call mma_deallocate(iAOtSO)
  nSOInf = 0
  nIrrep = 0

end subroutine SOAO_Info_Free

!***********************************************************************
!***********************************************************************

end module SOAO_Info
