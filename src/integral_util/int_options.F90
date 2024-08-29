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

module Int_Options

use Constants, only: Zero
use Definitions, only: wp, iwp

implicit none
private

integer(kind=iwp) :: iTOffs(8**3), Map4(4)
real(kind=wp) :: Disc = Zero, Disc_Mx = Zero, ExFac = Zero, Quad_ijkl = Zero, Thize = Zero
logical(kind=iwp) :: DoFock = .false., DoIntegrals = .true., FckNoClmb = .false., FckNoExch = .false., PreSch = .true., &
                     W2Disc = .false.

public :: Disc, Disc_Mx, DoFock, DoIntegrals, ExFac, FckNoClmb, FckNoExch, Init_Int_Options, iTOffs, Map4, PreSch, Quad_ijkl, &
          Thize, W2Disc

contains

subroutine Init_Int_Options()

  ! Default values
  DoFock = .false.
  DoIntegrals = .true.
  FckNoClmb = .false.
  FckNoExch = .false.
  PreSch = .true.
  W2Disc = .false.
  Disc = Zero
  Disc_Mx = Zero
  ExFac = Zero
  Quad_ijkl = Zero
  Thize = Zero

end subroutine Init_Int_Options

end module Int_Options
