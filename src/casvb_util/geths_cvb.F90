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

subroutine geths_cvb(arr,n)

use Definitions, only: iwp

implicit none
integer(kind=iwp), intent(inout) :: n
character(len=*), intent(inout) :: arr(n)
integer(kind=iwp) :: i, iaux(1), iret, j, lenarr

call geth_cvb(iaux,1)
n = iaux(1)
lenarr = len(arr(1))
do i=1,n
  do j=1,lenarr
    call geth_cvb(iaux,1)
    iret = iaux(1)
    arr(i)(j:j) = char(iret)
  end do
end do

return

end subroutine geths_cvb
