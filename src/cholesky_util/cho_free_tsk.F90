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

subroutine Cho_Free_Tsk(ID)

use Cholesky, only: Cho_Real_Par
use Definitions, only: iwp

implicit none
integer(kind=iwp) :: ID

if (Cho_Real_Par) then
  call Free_Tsk(ID)
else
  call Cho_Free_Tsk_(ID)
end if

end subroutine Cho_Free_Tsk
