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

use Definitions, only: wp, iwp
use Constants, only: Zero

implicit none
private

logical(kind=iwp) :: DoIntegrals = .true.
logical(kind=iwp) :: DoFock = .false.
logical(kind=iwp) :: FckNoClmb = .false.
logical(kind=iwp) :: FckNoExch = .false.
real(kind=wp) :: ExFac = Zero
real(kind=wp) :: Thize = Zero
real(kind=wp) :: Disc_Mx = Zero
real(kind=wp) :: Disc = Zero
logical(kind=iwp) :: W2Disc = .false.
logical(kind=iwp) :: PreSch = .true.
integer(kind=iwp) :: iTOffs(8**3)
real(kind=wp) :: Quad_ijkl = Zero
integer(kind=iwp) :: Map4(4)

public :: DoIntegrals, DoFock, FckNoClmb, FckNoExch, ExFac, Thize, W2Disc, PreSch, Disc_Mx, Disc, iTOffs, Quad_ijkl, Map4

end module Int_Options
