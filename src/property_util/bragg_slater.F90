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

function Bragg_Slater(iAtmNr)
! Bragg-Slater radii in Angstrom (J. C. Slater, Quantum Theory of
! Molecules and Solids, volume 2, table 3-1) up to Radon.

use Constants, only: Angstrom
use Definitions, only: wp, iwp, u6

implicit none
real(kind=wp) :: Bragg_Slater
integer(kind=iwp), intent(in) :: iAtmNr
real(kind=wp), parameter :: BS_Radii(102) = [ &
  0.25_wp,                                                0.25_wp,                                                         & ! 1-2
  1.45_wp,1.05_wp,0.85_wp,0.70_wp,0.65_wp,0.60_wp,0.50_wp,0.45_wp,                                                         & ! 3-10
  1.80_wp,1.50_wp,1.25_wp,1.10_wp,1.00_wp,1.00_wp,1.00_wp,1.00_wp,                                                         & ! 11-18
  2.20_wp,1.80_wp,                                                                                                         &
          1.60_wp,1.40_wp,1.35_wp,1.40_wp,1.40_wp,1.40_wp,1.35_wp,1.35_wp,1.35_wp,1.35_wp,                                 &
                  1.30_wp,1.25_wp,1.15_wp,1.15_wp,1.15_wp,1.15_wp,                                                         & ! 19-36
  2.35_wp,2.00_wp,                                                                                                         &
          1.80_wp,1.55_wp,1.45_wp,1.45_wp,1.35_wp,1.30_wp,1.35_wp,1.40_wp,1.60_wp,1.55_wp,                                 &
                  1.55_wp,1.45_wp,1.45_wp,1.40_wp,1.40_wp,1.40_wp,                                                         & ! 37-54
  2.60_wp,2.15_wp,                                                                                                         &
          1.95_wp,1.85_wp,1.85_wp,1.85_wp,1.85_wp,1.85_wp,1.85_wp,1.80_wp,1.75_wp,1.75_wp,1.75_wp,1.75_wp,1.75_wp,1.75_wp, &
          1.75_wp,1.55_wp,1.45_wp,1.35_wp,1.35_wp,1.30_wp,1.35_wp,1.35_wp,1.35_wp,1.50_wp,                                 &
                  1.90_wp,1.80_wp,1.60_wp,1.90_wp,1.90_wp,1.90_wp,                                                         & ! 55-86
  2.60_wp,2.15_wp,                                                                                                         &
          1.95_wp,1.80_wp,1.75_wp,1.75_wp,1.75_wp,1.75_wp,1.75_wp,1.75_wp,1.75_wp,1.75_wp,1.75_wp,1.75_wp,1.75_wp,1.75_wp  & ! 87-102
]

if (iAtmNr > size(BS_Radii)) then
  write(u6,*) 'Bragg-Slater: Too high atom number!'
  write(u6,*) 'iAtmNr=',iAtmNr
  call Quit_OnUserError()
end if
Bragg_Slater = BS_Radii(iAtmNr)/Angstrom

return

end function Bragg_Slater
