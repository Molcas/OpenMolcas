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

subroutine make_close_rvb()

use Definitions, only: iwp

implicit none
integer(kind=iwp) :: i, il, n
character(len=8) :: vec(11)
integer(kind=iwp), external :: find_lu

vec(1) = 'TMP01'
vec(2) = 'TMP02'
vec(3) = 'TMP03'
vec(4) = 'TMP04'
vec(5) = 'TMP05'
vec(6) = 'TMP06'
vec(7) = 'TMP07'
vec(8) = 'TMP08'
vec(9) = 'TMP09'
vec(10) = 'VBWFN'
il = 10
! Preassign some file names to identifiers :
do i=1,il
  n = find_lu(vec(i))
  if (n > 0) call daclos(n)
end do

end subroutine make_close_rvb
