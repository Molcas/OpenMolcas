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

subroutine findmx_cvb(vec,n,vmx,imx)

use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: n
real(kind=wp), intent(in) :: vec(n)
real(kind=wp), intent(out) :: vmx
integer(kind=iwp), intent(out) :: imx
integer(kind=iwp) :: i

if (n > 0) then
  imx = 1
  vmx = vec(1)
  do i=2,n
    if (vec(i) > vmx) then
      imx = i
      vmx = vec(i)
    end if
  end do
else
  imx = 0
  vmx = -1.0e20_wp
end if

return

end subroutine findmx_cvb
