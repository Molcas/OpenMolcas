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

subroutine RPA_Setup()

! Thomas Bondo Pedersen (CTCC,UiO), July 2013.
!
! Set up RPA calculation:
!    - process input
!    - read reference orbitals

implicit none
character(len=*), parameter :: SecNam = 'RPA_Setup'

! Define data in common blocks (dummy values).
call RPA_SetInc()

! Get and check integral representation from Runfile.
call RPA_SetIntegralRepresentation()
call RPA_CheckIntegralRepresentation()

! Pick up data from Runfile.
call RPA_RdRun()

! Process input.
call RPA_RdInp()

! Read orbitals and orbital energies.
call RPA_RdOrb()

! Postprocessing and print
call RPA_PPInp()

! Add info for testing
call RPA_Setup_Add_Info()

end subroutine RPA_Setup
