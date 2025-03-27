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
! Copyright (C) Anders Bernhardsson                                    *
!***********************************************************************

function Facti(R)

use Constants, only: One
use Definitions, only: wp

integer n, i, j, R
real*8 Facti

n = R
i = 1
if (n == 0) then
  Facti = One
  return
end if
do j=1,n
  i = i*j
end do
Facti = real(i,kind=wp)

return

end function Facti
