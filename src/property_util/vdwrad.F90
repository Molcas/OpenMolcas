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
! A function that returns the van der Waals radius of an element, when
! such a radius is available. The user should exercise some care and
! check if the radius is zero, in which case the radius has not been
! reported.
! References:
! Main group elements: J. Phys. Chem. A 113 (2009) 5806-5812 doi:10.1021/jp8111556
! Transition metals: Dalton Trans. 42 (2013) 8617-8636 doi:10.1039/c3dt50599e

use Constants, only: Angstrom
use Definitions, only: wp, iwp, u6

implicit none
real(kind=wp) :: vdwRad
integer(kind=iwp), intent(in) :: iAtmNr
real(kind=wp), parameter :: Radii(102) = [ &
  1.10_wp,                                                1.40_wp,                                                         & ! 1-2
  1.81_wp,1.53_wp,1.92_wp,1.70_wp,1.55_wp,1.52_wp,1.47_wp,1.54_wp,                                                         & ! 3-10
  2.27_wp,1.73_wp,1.84_wp,2.10_wp,1.80_wp,1.80_wp,1.75_wp,1.88_wp,                                                         & ! 11-18
  2.75_wp,2.31_wp,                                                                                                         &
          2.58_wp,2.46_wp,2.42_wp,2.45_wp,2.45_wp,2.44_wp,2.40_wp,2.40_wp,2.38_wp,2.39_wp,                                 &
                  1.87_wp,2.11_wp,1.85_wp,1.90_wp,1.83_wp,2.02_wp,                                                         & ! 19-36
  3.03_wp,2.49_wp,                                                                                                         &
          2.75_wp,2.52_wp,2.56_wp,2.45_wp,2.44_wp,2.46_wp,2.44_wp,2.15_wp,2.53_wp,2.49_wp,                                 &
                  1.93_wp,2.17_wp,2.06_wp,2.06_wp,1.98_wp,2.16_wp,                                                         & ! 37-54
  3.43_wp,2.68_wp,                                                                                                         &
          2.98_wp,2.88_wp,2.92_wp,2.95_wp,0.00_wp,2.90_wp,2.87_wp,2.83_wp,2.79_wp,2.87_wp,2.81_wp,2.83_wp,2.79_wp,2.80_wp, &
          2.74_wp,2.63_wp,2.53_wp,2.57_wp,2.49_wp,2.48_wp,2.41_wp,2.29_wp,2.32_wp,2.45_wp,                                 &
                  1.96_wp,2.02_wp,2.07_wp,1.97_wp,2.02_wp,2.20_wp,                                                         & ! 55-86
  3.48_wp,2.83_wp,                                                                                                         &
          2.80_wp,2.93_wp,2.88_wp,2.71_wp,2.82_wp,2.81_wp,2.83_wp,3.05_wp,3.40_wp,3.05_wp,2.70_wp,0.00_wp,0.00_wp,0.00_wp  & ! 87-102
]

if (iAtmNr > size(Radii)) then
  write(u6,*) 'vdWRad: Too high atom number!'
  write(u6,*) 'iAtmNr=',iAtmNr
  call Quit_OnUserError()
end if
vdWRad = Radii(iAtmNr)/Angstrom

return

end function vdWRad
