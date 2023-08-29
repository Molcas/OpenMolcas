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
! Copyright (C) 2005, Thomas Bondo Pedersen                            *
!***********************************************************************

subroutine Cho_InpMod(Mode)
!
! Thomas Bondo Pedersen, Jan. 2005.
!
! Purpose: modify Cholesky settings for
!          Mode = 'LOW ' : low-accuracy decomposition
!          Mode = 'MEDI' : medium-accuracy decomposition
!          Mode = 'HIGH' : high-accuracy decomposition
!          Mode = '1CCD' : set one-center approximation
!          (All other Mode-values are ignored.)

use Cholesky, only: Cho_1Center, ThrCom
use Definitions, only: wp

implicit none
character(len=4), intent(in) :: Mode
character(len=4) :: Mod2

Mod2 = Mode
call Upcase(Mod2)

if (Mod2(1:3) == 'LOW') then
  ThrCom = 1.0e-4_wp
else if (Mod2(1:4) == 'MEDI') then
  ThrCom = 1.0e-6_wp
else if (Mod2(1:4) == 'HIGH') then
  ThrCom = 1.0e-8_wp
else if (Mod2(1:4) == '1CCD') then
  Cho_1Center = .true.
end if

end subroutine Cho_InpMod
