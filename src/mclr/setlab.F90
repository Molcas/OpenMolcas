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
! Copyright (C) 1996, Anders Bernhardsson                              *
!***********************************************************************
      Subroutine SetLab(label,j)
      Implicit Real*8 (a-h,o-z)
      Character*(*) Label
      Logical no
      no=.true.
      do i=1,LEN(label)
      If (no.and.label(i:i).eq.' ') then
      Write(Label(i:i),'(I1)') j
      no=.false.
      end if
      end do
      Return
      end
