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

subroutine DefParo3v3Hlp2(i,Schem,Nomen)
! help routine to DefParo2v4, producing names of Disc files
! ex: Schem='XYZQ', i=1  ->  Nomen='XYZQ01'

use Definitions, only: iwp

implicit none
integer(kind=iwp), intent(in) :: i
character(len=4), intent(in) :: Schem
character(len=6), intent(out) :: Nomen

write(Nomen,'(a4,i2.2)') Schem,i

return

end subroutine DefParo3v3Hlp2
