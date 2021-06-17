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

implicit real*8(a-h,o-z)
parameter(nAtmNr=102)
real*8 Radii(nAtmNr)
save Radii
#include "angstr.fh"
data Radii/                                                &
  1.20d0,                                          1.40d0, & ! 1-2
  1.82d0,0.00d0,0.00d0,1.70d0,1.55d0,1.52d0,1.47d0,1.54d0, & ! 3-10
  2.27d0,1.73d0,0.00d0,2.10d0,1.80d0,1.80d0,1.75d0,1.88d0, & ! 10-18
  2.75d0,0.00d0,                                           &
         0.00d0,0.00d0,0.00d0,0.00d0,0.00d0,0.00d0,0.00d0,1.63d0,1.40d0,1.39d0, &
                1.87d0,0.00d0,1.85d0,1.90d0,1.85d0,2.02d0, & ! 19-36
  0.00d0,0.00d0,                                           &
         0.00d0,0.00d0,0.00d0,0.00d0,0.00d0,0.00d0,0.00d0,1.63d0,1.72d0,1.58d0, &
                1.93d0,2.17d0,0.00d0,2.06d0,1.98d0,2.16d0, & ! 37-54
  0.0d0,0.00d0,                                            &
         0.00d0,0.00d0,0.00d0,0.00d0,0.00d0,0.00d0,0.00d0,0.00d0,0.00d0,0.00d0,0.00d0,0.00d0,0.00d0,0.00d0, &
         0.00d0,0.00d0,0.00d0,0.00d0,0.00d0,0.00d0,0.00d0,1.75d0,1.66d0,1.55d0, &
                1.96d0,2.02d0,0.00d0,0.00d0,0.00d0,0.00d0, & ! 55-86
  0.0d0,0.00d0,                                            &
         0.00d0,0.00d0,0.00d0,1.86d0,0.00d0,0.00d0,0.00d0,0.00d0,0.00d0,0.00d0,0.00d0,0.00d0,0.00d0,0.00d0/ ! 87-102

Ang2au = 1.0d0/angstr
if (iAtmNr > nAtmNr) then
  write(6,*) 'vdWRad: Too high atom number!'
  write(6,*) 'iAtmNr=',iAtmNr
  call Quit_OnUserError()
end if
vdWRad = Radii(iAtmNr)*Ang2au

return

end function vdWRad
