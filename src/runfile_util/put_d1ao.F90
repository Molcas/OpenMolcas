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

subroutine Put_D1ao(D1ao,nD1ao)

use Definitions, only: wp, iwp

implicit none
integer(kind=iwp) :: nD1ao
real(kind=wp) :: D1ao(nD1ao)
character(len=24) :: Label

! Write the variational 1st order density matrix
! density matrix in AO/SO basis

Label = 'D1ao'
call Put_dArray(Label,D1ao,nD1ao)

return

end subroutine Put_D1ao
