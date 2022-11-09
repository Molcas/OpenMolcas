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
! Copyright (C) 1994, Per Ake Malmqvist                                *
!               2012, Thomas Bondo Pedersen                            *
!***********************************************************************

subroutine iTrnsps(n,m,a,b)
! Thomas Bondo Pedersen, August 2012: integer version of trnsps

use Definitions, only: iwp

implicit none
integer(kind=iwp), intent(in) :: n, m, a(n,m)
integer(kind=iwp), intent(out) :: b(m,n)
integer(kind=iwp) :: i

do i=1,n
  b(1:m,i) = a(i,1:m)
end do

return

end subroutine iTrnsps
