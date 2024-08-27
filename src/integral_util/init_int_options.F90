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

subroutine Init_Int_Options()

use Int_Options, only: DoIntegrals, DoFock, FckNoClmb, FckNoExch
use Int_Options, only: Thize, W2Disc, PreSch, Disc_Mx, Disc, Quad_ijkl
use Constants, only: Zero

implicit none

! Set variables in module Int_Options

DoIntegrals = .true.  ! Default value
DoFock = .false.      ! Default value
FckNoClmb = .false.   ! Default value
FckNoExch = .false.   ! Default value
Thize = Zero          ! Default value
W2Disc = .false.      ! Default value
PreSch = .true.       ! Default value
Disc_Mx = Zero        ! Default value
Disc = Zero           ! Default Value
Quad_ijkl = Zero      ! Default Value

end subroutine Init_Int_Options
