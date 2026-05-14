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

subroutine gethr_cvb(iarr,n)

use Definitions, only: iwp, RtoI

implicit none
integer(kind=iwp), intent(inout) :: n
integer(kind=iwp), intent(out) :: iarr(n)
integer(kind=iwp) :: iaux(1)

call geth_cvb(iaux,1)
n = iaux(1)
call geth_cvb(iarr,RtoI*n)

return

end subroutine gethr_cvb
