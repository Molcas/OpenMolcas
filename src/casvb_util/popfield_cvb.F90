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

subroutine popfield_cvb(ifc)
! IFC is field control:
! IFC=1 --> read to end of line only, no new line -- DISABLED in MOLCAS
! IFC=2 --> begin read from next line

use casvb_global, only: ifield, nfield, nfold
use Definitions, only: iwp

implicit none
integer(kind=iwp), intent(in) :: ifc
integer(kind=iwp) :: initpop = 0

if (initpop == 0) then
  ifield = 0
  nfield = 0
  nfold = 0
end if
initpop = 1

if ((ifield == nfield) .or. (ifc == 2)) then
  nfold = nfield
  call rdline_cvb(nfield)
  ifield = 1
else
  ! IFIELD > NFIELD will signify no read
  ifield = min(ifield+1,nfield+1)
end if

return

end subroutine popfield_cvb
