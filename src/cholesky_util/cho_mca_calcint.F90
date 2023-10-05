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

subroutine CHO_MCA_CALCINT(ISHLAB)
!
! Purpose: calculate qualified integral columns from
!          shell pair distribution (**|ISHLA ISHLB).

use Cholesky, only: IFCSEW
use Definitions, only: iwp

implicit none
integer(kind=iwp), intent(in) :: ISHLAB
character(len=*), parameter :: SECNAM = 'CHO_MCA_CALCINT'

if (IFCSEW == 1) then ! store full shell quadruple
  call CHO_MCA_CALCINT_1(ISHLAB)
else if (IFCSEW == 2) then ! store only reduced sets
  call CHO_MCA_CALCINT_2(ISHLAB)
else ! this is an error
  call CHO_QUIT('IFCSEW out of bounds in '//SECNAM,105)
end if

end subroutine CHO_MCA_CALCINT
