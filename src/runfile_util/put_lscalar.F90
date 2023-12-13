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

subroutine Put_lScalar(Label,Logc)

use Definitions, only: iwp

implicit none
character(len=*), intent(in) :: Label
logical(kind=iwp), intent(in) :: Logc

call Put_iScalar(Label,merge(1,0,Logc))

return

end subroutine Put_lScalar
