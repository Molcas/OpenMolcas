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

implicit real*8(a-h,o-z)
parameter(nAtmNr=102)
real*8 BS_Radii(nAtmNr)
save BS_Radii
#include "angstr.fh"
data BS_Radii/                                                     &
  0.25d0,0.25d0,                                                   &
  1.45d0,1.05d0,0.85d0,0.70d0,0.65d0,0.60d0,0.50d0,0.45d0,         &
  1.80d0,1.50d0,1.25d0,1.10d0,1.00d0,1.00d0,1.00d0,1.00d0,         &
  2.20d0,1.80d0,                                                   & ! 19-36
         1.60d0,1.40d0,1.35d0,1.40d0,1.40d0,                       &
                         1.40d0,1.35d0,1.35d0,1.35d0,1.35d0,       &
                         1.30d0,1.25d0,1.15d0,1.15d0,1.15d0,1.15d0,&
  2.35d0,2.00d0,                                                   & ! 37-54
         1.80d0,1.55d0,1.45d0,1.45d0,1.35d0,                       &
                         1.30d0,1.35d0,1.40d0,1.60d0,1.55d0,       &
                         1.55d0,1.45d0,1.45d0,1.40d0,1.40d0,1.40d0,&
  2.6d0,2.15d0,                                                    & ! 55-86
         1.95d0,1.85d0,1.85d0,1.85d0,1.85d0,1.85d0,1.85d0,         &
                1.80d0,1.75d0,1.75d0,1.75d0,1.75d0,1.75d0,1.75d0,  &
         1.75d0,1.55d0,1.45d0,1.35d0,1.35d0,                       &
                         1.30d0,1.35d0,1.35d0,1.35d0,1.50d0,       &
                         1.90d0,1.80d0,1.60d0,1.90d0,1.90d0,1.90d0,&
  2.6d0,2.15d0,                                                    & ! 87-102
         1.95d0,1.80d0,1.75d0,1.75d0,1.75d0,1.75d0,1.75d0,         &
                1.75d0,1.75d0,1.75d0,1.75d0,1.75d0,1.75d0,1.75d0/

Ang2au = 1.0d0/angstr
if (iAtmNr > nAtmNr) then
  write(6,*) 'Bragg-Slater: Too high atom number!'
  write(6,*) 'iAtmNr=',iAtmNr
  call Quit_OnUserError()
end if
Bragg_Slater = BS_Radii(iAtmNr)*Ang2au

return

end function Bragg_Slater
