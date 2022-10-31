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
integer(kind=iwp) :: n, m, a(n,m), b(m,n)
integer(kind=iwp) :: i, j

do i=1,n
  do j=1,m
    b(j,i) = a(i,j)
  end do
end do

return

end subroutine iTrnsps
