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
! Copyright (C) 2013, Thomas Bondo Pedersen                            *
!***********************************************************************

subroutine RPA_PPInp()

! Thomas Bondo Pedersen (CTCC,UiO), July 2013.
!
! Input postprocessing.

use RPA_globals, only: dRPA, Reference, RPAModel, SOSEX

implicit none
character(len=*), parameter :: SecNam = 'RPA_PPInp'

! set RPAModel
if (dRPA) then
  if (SOSEX) then
    RPAModel = 'SOSX@'//Reference(1:3)
  else
    RPAModel = 'dRPA@'//Reference(1:3)
  end if
else
  ! this should never happen
  call RPA_Warn(3,SecNam//': internal error [RPAModel]')
  RPAModel = 'None@Non'
end if

! freeze orbitals
call RPA_Freezer()

! print config after input processing
call RPA_PrInp()

end subroutine RPA_PPInp
