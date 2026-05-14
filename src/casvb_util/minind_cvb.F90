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

function minind_cvb(iminor,i_dim,nel,ixmin)

use Definitions, only: iwp

implicit none
integer(kind=iwp) :: minind_cvb
integer(kind=iwp), intent(in) :: i_dim, iminor(i_dim), nel, ixmin(0:nel,0:i_dim)
integer(kind=iwp) :: i

minind_cvb = 1
do i=1,i_dim
  minind_cvb = minind_cvb+ixmin(iminor(i)-1,i)
end do

return

end function minind_cvb
