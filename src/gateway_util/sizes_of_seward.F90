!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!***********************************************************************

module Sizes_of_Seward

use define_af, only: iTabMx
use Definitions, only: iwp

implicit none
private

public :: S, Size_Dmp, Size_Get

integer(kind=iwp), parameter :: nLen = 16+2*iTabMx ! number of elements
type Sizes_of_Stuff
  integer(kind=iwp) :: m2Max = 0
  integer(kind=iwp) :: mCentr = 0
  integer(kind=iwp) :: mCentr_Aux = 0
  integer(kind=iwp) :: mCentr_Frag = 0
  integer(kind=iwp) :: Mx_mdc = 0
  integer(kind=iwp) :: Mx_Shll = 0
  integer(kind=iwp) :: n2Tot = 0
  integer(kind=iwp) :: jMax = 5
  integer(kind=iwp) :: MaxPrm(0:iTabMx) = 0
  integer(kind=iwp) :: MaxBas(0:iTabMx) = 0
  integer(kind=iwp) :: nDim = 0
  integer(kind=iwp) :: nShlls = 0
  integer(kind=iwp) :: Max_Center = 15
  integer(kind=iwp) :: kCentr = 0
  integer(kind=iwp) :: nMltpl = 2
  integer(kind=iwp) :: iAngMx = -1
end type Sizes_of_Stuff

type(Sizes_of_Stuff) :: S

contains

subroutine Size_Dmp()

  use stdalloc, only: mma_allocate, mma_deallocate

  integer(kind=iwp) :: i
  integer(kind=iwp), allocatable :: iDmp(:)

  call mma_allocate(iDmp,nLen,Label='iDmp')

  i = 0
  iDmp(i+1) = S%m2Max
  i = i+1
  iDmp(i+1) = S%mCentr
  i = i+1
  iDmp(i+1) = S%mCentr_Aux
  i = i+1
  iDmp(i+1) = S%mCentr_Frag
  i = i+1
  iDmp(i+1) = S%Mx_mdc
  i = i+1
  iDmp(i+1) = S%Mx_Shll
  i = i+1
  iDmp(i+1) = S%n2Tot
  i = i+1
  iDmp(i+1) = S%jMax
  i = i+1
  iDmp(i+1:i+iTabMx+1) = S%MaxPrm(0:iTabMx)
  i = i+iTabMx+1
  iDmp(i+1:i+iTabMx+1) = S%MaxBas(0:iTabMx)
  i = i+iTabMx+1
  iDmp(i+1) = S%nDim
  i = i+1
  iDmp(i+1) = S%nShlls
  i = i+1
  iDmp(i+1) = S%Max_Center
  i = i+1
  iDmp(i+1) = S%kCentr
  i = i+1
  iDmp(i+1) = S%nMltpl
  i = i+1
  iDmp(i+1) = S%iAngMx
  i = i+1

  call Put_iArray('Sizes',iDmp,nLen)
  call mma_deallocate(iDmp)

end subroutine Size_Dmp

subroutine Size_Get()

  use stdalloc, only: mma_allocate, mma_deallocate
  use Definitions, only: u6

  integer(kind=iwp) :: i, Len2
  logical(kind=iwp) :: Found
  integer(kind=iwp), allocatable :: iDmp(:)

  call mma_allocate(iDmp,nLen,Label='iDmp')

  call Qpg_iArray('Sizes',Found,Len2)
  if (.not. Found) then
    write(u6,*) 'Size_Get: Sizes not found.'
    call Abend()
  end if
  if (nLen /= Len2) then
    write(u6,*) 'Size_Get: nLen /= Len2.'
    call Abend()
  end if
  call Get_iArray('Sizes',iDmp,nLen)

  i = 0
  S%m2Max = iDmp(i+1)
  i = i+1
  S%mCentr = iDmp(i+1)
  i = i+1
  S%mCentr_Aux = iDmp(i+1)
  i = i+1
  S%mCentr_Frag = iDmp(i+1)
  i = i+1
  S%Mx_mdc = iDmp(i+1)
  i = i+1
  S%Mx_Shll = iDmp(i+1)
  i = i+1
  S%n2Tot = iDmp(i+1)
  i = i+1
  S%jMax = iDmp(i+1)
  i = i+1
  S%MaxPrm(0:iTabMx) = iDmp(i+1:i+1+iTabMx)
  i = i+iTabMx+1
  S%MaxBas(0:iTabMx) = iDmp(i+1:i+1+iTabMx)
  i = i+iTabMx+1
  S%nDim = iDmp(i+1)
  i = i+1
  S%nShlls = iDmp(i+1)
  i = i+1
  S%Max_Center = iDmp(i+1)
  i = i+1
  S%kCentr = iDmp(i+1)
  i = i+1
  S%nMltpl = iDmp(i+1)
  i = i+1
  S%iAngMx = iDmp(i+1)
  i = i+1

  call mma_deallocate(iDmp)

end subroutine Size_Get

end module Sizes_of_Seward
