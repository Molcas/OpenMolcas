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
! Copyright (C) 1996-2006, Thorstein Thorsteinsson                     *
!               1996-2006, David L. Cooper                             *
!***********************************************************************
!**************************************
!** Low-level input parsing routines **
!**************************************

subroutine pushfield_cvb()

use casvb_global, only: ifield, nfield, nfold

implicit none

if ((ifield == 1) .or. (nfield == -1)) then
  call pushline_cvb()
  ifield = nfold
  nfield = nfold
else
  ifield = ifield-1
end if

return

end subroutine pushfield_cvb
