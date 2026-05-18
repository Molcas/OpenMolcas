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

module fake_ga

use stdalloc, only: mma_allocate, mma_deallocate
use constants, only: Zero
use definitions, only: iwp, wp, u6

implicit none
private

integer(kind=iwp), parameter :: max_ga_arrays = 10
real(kind=wp) DBL_MB(2)
integer(kind=iwp) :: iga_arrays = 0
type ga_type
  integer(kind=iwp) :: g_a           ! global array handler
  real(kind=wp), allocatable :: A(:) ! local array
  integer(kind=iwp) :: index         ! Array(1)=DBL_MB(index)
  integer(kind=iwp) :: iLow = 1
  integer(kind=iwp) :: iHi
  integer(kind=iwp) :: jLow = 1
  integer(kind=iwp) :: jHi
  integer(kind=iwp) :: Length
end type ga_type

type(ga_type) :: GA_arrays(max_ga_arrays)

public :: GA_arrays, Allocate_GA_Array, Deallocate_GA_Array, DBL_MB

contains

integer(kind=iwp) function Allocate_GA_Array(nSize,Label) result(lg_A)

  integer(kind=iwp), intent(in) :: nSize
  character(len=*), intent(in) :: Label
  integer(kind=iwp) i

  lg_A = 0
  do i=1,max_ga_arrays
    if (.not. allocated(GA_arrays(i)%A)) then
      iga_arrays = iga_arrays+1
      call mma_allocate(GA_arrays(i)%A,nSize,Label=Label)
      lg_A = i
      GA_Arrays(lg_A)%A(:) = Zero
      return
    end if
  end do
  write(u6,*) 'To many GA_arrys, increase max_ga_arrays.'
  call abend()

end function Allocate_GA_Array

subroutine Deallocate_GA_Array(lg_A)

  integer(kind=iwp), intent(inout) :: lg_A

  call mma_deallocate(GA_Arrays(lg_A)%A)
  iga_arrays = iga_arrays-1
  lg_A = 0

end subroutine Deallocate_GA_Array

end module fake_ga
