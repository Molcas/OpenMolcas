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

function JSUNP(INTSYM,L)

use Definitions, only: iwp

implicit none
integer(kind=iwp) :: JSUNP
integer(kind=iwp), intent(in) :: INTSYM(*), L
integer(kind=iwp) :: I, J

! L = (I-1)*10 + J+1  [ J = 0-9 ]
I = (L+9)/10
J = mod(L-1,10)
JSUNP = 1+ibits(INTSYM(I),27-3*J,3)

return

end function JSUNP
