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

!IFG trivial
subroutine cblank_cvb(a,ndim)

use Definitions, only: iwp

implicit none
character(len=*) :: a
integer(kind=iwp) :: ndim
integer(kind=iwp) :: i
character, parameter :: blank = ' '

do i=1,ndim
  a(i:i) = blank
end do

return

end subroutine cblank_cvb
