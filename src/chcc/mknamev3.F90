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

subroutine MkNameV3(i,j,k,Schem,Nomen)
! help routine to ReaW3, producing name of V3file
! ex: Schem='XY', i=1, j=3, k=5 ->  Nomen='XY010305'

use Definitions, only: iwp

implicit none
integer(kind=iwp), intent(in) :: i, j, k
character(len=2), intent(in) :: Schem
character(len=8), intent(out) :: Nomen

write(Nomen,'(a2,3(i2.2))') Schem,i,j,k

return

end subroutine MkNameV3
