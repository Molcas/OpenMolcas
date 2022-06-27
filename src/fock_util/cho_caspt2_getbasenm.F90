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

subroutine Cho_caspt2_GetBaseNm(BaseNm,iTyp)

use Definitions, only: iwp

implicit none
character(len=3), intent(out) :: BaseNm
integer(kind=iwp), intent(in) :: iTyp

if (iTyp == 1) then
  BaseNm = '_PI'
else if (iTyp == 2) then
  BaseNm = '_PW'
else if (iTyp == 3) then
  BaseNm = '_CD'
else
  BaseNm = '_un'
end if

end subroutine Cho_caspt2_GetBaseNm
