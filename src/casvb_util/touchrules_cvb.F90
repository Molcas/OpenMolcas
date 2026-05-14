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
! Copyright (C) 1996-2006, Thorstein Thorsteinsson                     *
!               1996-2006, David L. Cooper                             *
!***********************************************************************

subroutine touchrules_cvb(chr)

implicit none
character(len=*), intent(in) :: chr

select case (chr)
  case ('CI-ORBS')
    call clearcnt_cvb(1)
  case ('CI-CVB')
    call clearcnt_cvb(2)
  case ('CI-ALL')
    call clearcnt_cvb(3)
end select

return

end subroutine touchrules_cvb
