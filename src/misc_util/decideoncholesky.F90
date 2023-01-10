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

subroutine DecideonCholesky(DoCholesky)

use Definitions, only: iwp

implicit none
logical(kind=iwp), intent(out) :: DoCholesky
integer(kind=iwp) :: iOption

call Get_iScalar('System BitSwitch',iOption)

! this is a logical variable (logical function)
DoCholesky = btest(iOption,9)

return

end subroutine DecideonCholesky
