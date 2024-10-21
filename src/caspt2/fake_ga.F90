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

Real*8 DBL_MB(2)

Type ga_type
     Integer :: g_a           ! global array handler
     Real*8, allocatable  :: Array(:)      ! local array
     Integer :: index         ! Array(1)=DBL_MB(index)
     Integer :: iLow = 1
     Integer :: iHi
     Integer :: jLow = 1
     Integer :: jHi
     Integer :: Length
End Type ga_type

Integer, parameter :: max_ga_arrays=10
Integer :: iga_arrays_=0

Type (ga_type) :: GA_arrays(max_ga_arrays)

End Module fake_ga
