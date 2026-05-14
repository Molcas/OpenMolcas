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

subroutine DefParo3v3Hlp1(i,j,Schem,Nomen)
! help routine to DefParo2v4, producing names of Disc files
! ex: Schem='XY', i=1, j=3  ->  Nomen='XY0103'
! N.B. suspendovana rutina

use Definitions, only: iwp

implicit none
integer(kind=iwp), intent(in) :: i, j
character(len=2), intent(in) :: Schem
character(len=6), intent(out) :: Nomen

write(Nomen,'(a2,2(i2.2))') Schem,i,j

return

end subroutine DefParo3v3Hlp1
