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

subroutine seths_cvb(arr,n)

use Definitions, only: iwp

implicit none
integer(kind=iwp), intent(in) :: n
character(len=*), intent(in) :: arr(n)
integer(kind=iwp) :: i, j, lenarr

call seth_cvb([n],1)
lenarr = len(arr)
do i=1,n
  do j=1,lenarr
    call seth_cvb([ichar(arr(i)(j:j))],1)
  end do
end do

return

end subroutine seths_cvb
