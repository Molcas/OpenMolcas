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
! Copyright (C) 2004, Thomas Bondo Pedersen                            *
!***********************************************************************

subroutine ChoMP2_GetBaseNm(BaseNm,iTyp)

use Definitions, only: iwp

implicit none
character(len=3), intent(out) :: BaseNm
integer(kind=iwp), intent(in) :: iTyp

if (iTyp == 1) then
  BaseNm = '_AI'
else if (iTyp == 2) then
  BaseNm = '_CD'
else
  BaseNm = '_un'
end if

end subroutine ChoMP2_GetBaseNm
