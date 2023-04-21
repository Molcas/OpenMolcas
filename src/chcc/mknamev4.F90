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

subroutine MkNameV4(i,j,k,l,Schem,Nomen)
! help routine to ReaW4, producing name of V4file
! ex: Schem='XY', i=1, j=3, k=5, l=07 ->  Nomen='XY01030507'

use Definitions, only: iwp

implicit none
integer(kind=iwp), intent(in) :: i, j, k, l
character(len=2), intent(in) :: Schem
character(len=10), intent(out) :: Nomen

write(Nomen,'(a2,4(i2.2))') Schem,i,j,k,l

return

end subroutine MkNameV4
