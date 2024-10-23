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
Module fake_ga
use stdalloc, only: mma_allocate, mma_deallocate
Private

Real*8 DBL_MB(2)

Type ga_type
     Integer :: g_a           ! global array handler
     Real*8, allocatable  :: A(:)      ! local array
     Integer :: index         ! Array(1)=DBL_MB(index)
     Integer :: iLow = 1
     Integer :: iHi
     Integer :: jLow = 1
     Integer :: jHi
     Integer :: Length
End Type ga_type

Integer, parameter :: max_ga_arrays=10
Integer :: iga_arrays=0

Type (ga_type) :: GA_arrays(max_ga_arrays)

Public :: GA_arrays, Allocate_GA_Array, Deallocate_GA_Array, DBL_MB

Contains
Integer Function Allocate_GA_Array(nSize,Label) result(lg_A)
Implicit None
Integer, Intent(In):: nSize
Character(LEN=*), Intent(In):: Label

Integer i

lg_A=0
Do i = 1, max_ga_arrays
   If (.Not.Allocated(GA_arrays(i)%A)) Then
      iga_arrays=iga_arrays+1
      Call mma_allocate(GA_arrays(i)%A,nSize,Label=Label)
      lg_A=i
      GA_Arrays(lg_A)%A(:)=0.0D0
      Return
   End If
End Do
Write (6,*) 'To many GA_arrys, increase max_ga_arrays.'
Call abend()
End Function Allocate_GA_Array

Subroutine Deallocate_GA_Array(lg_A)
Integer, Intent(InOut):: lg_A
Call mma_deallocate(GA_Arrays(lg_A)%A)
iga_arrays=iga_arrays-1
lg_A=0
End Subroutine Deallocate_GA_Array
End Module fake_ga
