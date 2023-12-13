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

subroutine rea1(lun,length,A)
! nacitane bloku dat z Sekvencneho suboru, ak koniec tak
! rewind a citanie znova od zaciatku
! Toto je ozaj ditry

use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: lun, length
real(kind=wp), intent(out) :: A(length)
integer(kind=iwp) :: istatus

read(lun,iostat=istatus) A
if (istatus < 0) then
  rewind(lun)
  read(lun) A
end if

return

end subroutine rea1
