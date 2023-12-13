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

function vdWRad(iAtmNr)
! A function that returns the van der Waals radie of an element, when
! such a radie is available. The user should exercise some care and
! check if the radie is zero, in which case the radie has not been
! reported. Reference: Bondi, J.Phys.Chem. 68 (1964) 441.

use Constants, only: Angstrom
use Definitions, only: wp, iwp, u6

implicit none
real(kind=wp) :: vdwRad
integer(kind=iwp), intent(in) :: iAtmNr
real(kind=wp), parameter :: Radii(102) = [ &
  1.20_wp,                                                1.40_wp,                                                         & ! 1-2
  1.82_wp,0.00_wp,0.00_wp,1.70_wp,1.55_wp,1.52_wp,1.47_wp,1.54_wp,                                                         & ! 3-10
  2.27_wp,1.73_wp,0.00_wp,2.10_wp,1.80_wp,1.80_wp,1.75_wp,1.88_wp,                                                         & ! 11-18
  2.75_wp,0.00_wp,                                                                                                         &
          0.00_wp,0.00_wp,0.00_wp,0.00_wp,0.00_wp,0.00_wp,0.00_wp,1.63_wp,1.40_wp,1.39_wp,                                 &
                  1.87_wp,0.00_wp,1.85_wp,1.90_wp,1.85_wp,2.02_wp,                                                         & ! 19-36
  0.00_wp,0.00_wp,                                                                                                         &
          0.00_wp,0.00_wp,0.00_wp,0.00_wp,0.00_wp,0.00_wp,0.00_wp,1.63_wp,1.72_wp,1.58_wp,                                 &
                  1.93_wp,2.17_wp,0.00_wp,2.06_wp,1.98_wp,2.16_wp,                                                         & ! 37-54
  0.00_wp,0.00_wp,                                                                                                         &
          0.00_wp,0.00_wp,0.00_wp,0.00_wp,0.00_wp,0.00_wp,0.00_wp,0.00_wp,0.00_wp,0.00_wp,0.00_wp,0.00_wp,0.00_wp,0.00_wp, &
          0.00_wp,0.00_wp,0.00_wp,0.00_wp,0.00_wp,0.00_wp,0.00_wp,1.75_wp,1.66_wp,1.55_wp,                                 &
                  1.96_wp,2.02_wp,0.00_wp,0.00_wp,0.00_wp,0.00_wp,                                                         & ! 54-86
  0.00_wp,0.00_wp,                                                                                                         &
          0.00_wp,0.00_wp,0.00_wp,1.86_wp,0.00_wp,0.00_wp,0.00_wp,0.00_wp,0.00_wp,0.00_wp,0.00_wp,0.00_wp,0.00_wp,0.00_wp  & ! 87-102
]

if (iAtmNr > size(Radii)) then
  write(u6,*) 'vdWRad: Too high atom number!'
  write(u6,*) 'iAtmNr=',iAtmNr
  call Quit_OnUserError()
end if
vdWRad = Radii(iAtmNr)/Angstrom

return

end function vdWRad
