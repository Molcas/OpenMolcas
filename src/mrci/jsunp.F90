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
integer(kind=iwp) :: INTSYM(*), L

JSUNP = 1+mod(INTSYM((L+9)/10)/(2**(27-3*mod(L-1,10))),8)

return

end function JSUNP
