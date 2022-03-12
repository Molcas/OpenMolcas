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
!-------------------------------------------------------------------------*
! A subroutine that emulates the len_trim of later Fortran versions, but  *
! that is missing in some Fortran 77 compilers.                           *
!-------------------------------------------------------------------------*
      Integer Function Len_TrimAO(String)
      Character*(*) String
      Do 15,i=Len(String),1,-1
        If(String(i:i).ne.' ') Go To 20
15    Continue
20    Continue
      Len_TrimAO=i
      Return
      End
