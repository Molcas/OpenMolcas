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
<<<<<<< HEAD
Module fake_ga
use stdalloc, only: mma_allocate, mma_deallocate
use constants, only:Zero
use definitions, only: iwp, wp, u6
Private

Real(kind=wp) DBL_MB(2)

Type ga_type
     Integer(kind=iwp) :: g_a           ! global array handler
     Real(kind=wp), allocatable  :: A(:)      ! local array
     Integer(kind=iwp) :: index         ! Array(1)=DBL_MB(index)
     Integer(kind=iwp) :: iLow = 1
     Integer(kind=iwp) :: iHi
     Integer(kind=iwp) :: jLow = 1
     Integer(kind=iwp) :: jHi
     Integer(kind=iwp) :: Length
End Type ga_type

Integer(kind=iwp), parameter :: max_ga_arrays=10
Integer(kind=iwp) :: iga_arrays=0

Type (ga_type) :: GA_arrays(max_ga_arrays)

Public :: GA_arrays, Allocate_GA_Array, Deallocate_GA_Array, DBL_MB

Contains
Integer Function Allocate_GA_Array(nSize,Label) result(lg_A)
Implicit None
Integer(kind=iwp), Intent(In):: nSize
Character(LEN=*), Intent(In):: Label

Integer(kind=iwp) i

lg_A=0
Do i = 1, max_ga_arrays
   If (.Not.Allocated(GA_arrays(i)%A)) Then
      iga_arrays=iga_arrays+1
      Call mma_allocate(GA_arrays(i)%A,nSize,Label=Label)
      lg_A=i
      GA_Arrays(lg_A)%A(:)=Zero
      Return
   End If
End Do
Write (u6,*) 'To many GA_arrys, increase max_ga_arrays.'
Call abend()
End Function Allocate_GA_Array

Subroutine Deallocate_GA_Array(lg_A)
Integer(kind=iwp), Intent(InOut):: lg_A
Call mma_deallocate(GA_Arrays(lg_A)%A)
iga_arrays=iga_arrays-1
lg_A=0
End Subroutine Deallocate_GA_Array
End Module fake_ga
=======

module fake_ga

use Data_Structures, only: Alloc1DArray_Type
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero
use Definitions, only: wp, iwp, u6

implicit none
private

integer(kind=iwp), parameter :: max_ga_arrays = 10

integer(kind=iwp) :: iga_arrays = 0
real(kind=wp) :: DBL_MB(2)
type(Alloc1DArray_Type) :: GA_arrays(max_ga_arrays)

public :: GA_arrays, Allocate_GA_Array, Deallocate_GA_Array, DBL_MB

contains

function Allocate_GA_Array(nSize,Label) result(lg_A)

  integer(kind=iwp) :: lg_A
  integer(kind=iwp), intent(in) :: nSize
  character(len=*), intent(in) :: Label
  integer(kind=iwp) :: i

  lg_A = 0
  do i=1,max_ga_arrays
    if (.not. allocated(GA_arrays(i)%A)) exit
  end do
  if (i <= max_ga_arrays) then
    iga_arrays = iga_arrays+1
    call mma_allocate(GA_arrays(i)%A,nSize,Label=Label)
    lg_A = i
    GA_Arrays(lg_A)%A(:) = Zero
  else
    write(u6,*) 'To many GA_arrays, increase max_ga_arrays.'
    call abend()
  end if

end function Allocate_GA_Array

subroutine Deallocate_GA_Array(lg_A)

  integer(kind=iwp), intent(inout) :: lg_A

  call mma_deallocate(GA_Arrays(lg_A)%A)
  iga_arrays = iga_arrays-1
  lg_A = 0

end subroutine Deallocate_GA_Array

end module fake_ga
>>>>>>> upstream-openmolcas/master
