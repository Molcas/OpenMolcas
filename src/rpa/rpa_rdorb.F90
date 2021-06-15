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

subroutine RPA_RdOrb()

! Thomas Bondo Pedersen (CTCC,UiO), July 2013.
!
! Read orbitals and orbital energies from InpOrb or from Runfile.

use RPA_globals, only: LumOrb

implicit none
character(len=*), parameter :: SecNam = 'RPA_RdOrb'

if (LumOrb) then
  ! read from InpOrb
  call RPA_RdOrb_FromInpOrb()
else
  ! read from Runfile
  call RPA_RdOrb_FromRunfile()
end if

end subroutine RPA_RdOrb
