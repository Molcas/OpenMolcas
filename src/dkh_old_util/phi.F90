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
! Copyright (C) 1984,1986, Bernd Artur Hess                            *
!***********************************************************************

function PHI(M,N)

use crelop, only: GA
use Constants, only: Zero, Two
use Definitions, only: wp, iwp

implicit none
real(kind=wp) :: PHI
integer(kind=iwp), intent(in) :: M, N

if ((mod(N,2) == 1) .or. (mod(M,2) == 1)) then
  PHI = Zero
else
  PHI = Two*GA(M+1)*GA(N+1)/GA(M+N+2)
end if

return

end function PHI
