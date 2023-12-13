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
      Subroutine xRlsMem_Ints()
      use k2_arrays, only: XMem
      implicit none
      If (.NOT.XMem) Then
!        Write (6,*)
!        Write (6,*) 'xRlsMem_Ints:',
!    &               'No external scratch handling to deactivate!'
!        Write (6,*)
      Else
!        Write (6,*)
!        Write (6,*) 'xRlsMem_Ints:','External allocation deactivate!'
!        Write (6,*)
         XMem=.False.
         Call RlsMem_Ints()
!        Write (6,*) 'xRlsMem_Ints:','External allocation released!'
      End If
!
      End Subroutine xRlsMem_Ints
