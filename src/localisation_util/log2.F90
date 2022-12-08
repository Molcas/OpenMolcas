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
! Copyright (C) Francesco Aquilante                                    *
!***********************************************************************

function Log2(n)

use Definitions, only: iwp

implicit none
integer(kind=iwp) :: Log2
integer(kind=iwp), intent(in) :: n
integer(kind=iwp) :: m

m = n
Log2 = 0
do while (m > 1)
  m = m/2
  Log2 = Log2+1
end do

return

end function Log2
