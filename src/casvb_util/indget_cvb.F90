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

function indget_cvb(iminor,i_dim,nel,ixmin)

use Definitions, only: iwp

implicit none
integer(kind=iwp) :: indget_cvb
integer(kind=iwp), intent(in) :: nel, iminor(nel), i_dim, ixmin(0:nel,0:i_dim)
integer(kind=iwp) :: i, iacc

indget_cvb = 1
iacc = 0
do i=1,nel
  if (iminor(i) == 1) then
    iacc = iacc+1
    indget_cvb = indget_cvb+ixmin(i-1,iacc)
  end if
end do

return

end function indget_cvb
