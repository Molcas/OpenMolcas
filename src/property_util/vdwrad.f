************************************************************************
* This file is part of OpenMolcas.                                     *
*                                                                      *
* OpenMolcas is free software; you can redistribute it and/or modify   *
* it under the terms of the GNU Lesser General Public License, v. 2.1. *
* OpenMolcas is distributed in the hope that it will be useful, but it *
* is provided "as is" and without any express or implied warranties.   *
* For more details see the full text of the license in the file        *
* LICENSE or in <http://www.gnu.org/licenses/>.                        *
************************************************************************
      Function vdWRad(iAtmNr)
*
*
*-- A function that returns the van der Waals radie of an element, when
*   such a radie is available. The user should exercise some care and
*   check if the radie is zero, in which case the radie has not been
*   reported. Reference: Bondi, J.Phys.Chem. 68 (1964) 441.
*
      Implicit Real*8 (a-h,o-z)
      Parameter (nAtmNr=102)
      Real*8 Radii(nAtmNr)
      save Radii
#include "angstr.fh"
      Data Radii/
     & 1.20D0,                                                   1.40D0,
     & 1.82D0,0.00D0,         0.00D0,1.70D0,1.55D0,1.52D0,1.47D0,1.54D0,
     & 2.27D0,1.73D0,         0.00D0,2.10D0,1.80D0,1.80D0,1.75D0,1.88D0,
     & 2.75D0,0.00D0,                                           ! 19-36
     &        0.00D0,0.00D0,0.00D0,0.00D0,0.00D0,
     &                        0.00D0,0.00D0,1.63D0,1.40D0,1.39D0,
     &                        1.87D0,0.00D0,1.85D0,1.90D0,1.85D0,2.02D0,
     & 0.00D0,0.00D0,                                           ! 37-54
     &        0.00D0,0.00D0,0.00D0,0.00D0,0.00D0,
     &                        0.00D0,0.00D0,1.63D0,1.72D0,1.58D0,
     &                        1.93D0,2.17D0,0.00D0,2.06D0,1.98D0,2.16D0,
     & 0.0D0,0.00D0,                                           ! 55-86
     &        0.00D0,0.00D0,0.00D0,0.00D0,0.00D0,0.00D0,0.00D0,
     &               0.00D0,0.00D0,0.00D0,0.00D0,0.00D0,0.00D0,0.00D0,
     &        0.00D0,0.00D0,0.00D0,0.00D0,0.00D0,
     &                        0.00D0,0.00D0,1.75D0,1.66D0,1.55D0,
     &                        1.96D0,2.02D0,0.00D0,0.00D0,0.00D0,0.00D0,
     & 0.0D0,0.00D0,                                          ! 87-102
     &        0.00D0,0.00D0,0.00D0,1.86D0,0.00D0,0.00D0,0.00D0,
     &               0.00D0,0.00D0,0.00D0,0.00D0,0.00D0,0.00D0,0.00D0/
*
      Ang2au=1.0D0/angstr
      If (iAtmNr.gt.nAtmNr) Then
        Write(6,*) 'vdWRad: Too high atom number!'
        Write(6,*) 'iAtmNr=',iAtmNr
        Call Quit_OnUserError
      End If
      vdWRad=Radii(iAtmNr)*Ang2au
      Return
      End
