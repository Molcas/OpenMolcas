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
      Function Bragg_Slater(iAtmNr)
      Implicit Real*8 (a-h,o-z)
*     Bragg-Slater radii in Angstrom (J. C. Slater, Quantum Theory of
*     Molecules and Solids, volume 2, table 3-1) up to Radon.
      Parameter (nAtmNr=102)
      Real*8 BS_Radii(nAtmNr)
      save BS_Radii
#include "angstr.fh"
      Data BS_Radii/
     & 0.25D0,                                                   0.25D0,
     & 1.45D0,1.05D0,         0.85D0,0.70D0,0.65D0,0.60D0,0.50D0,0.45D0,
     & 1.80D0,1.50D0,         1.25D0,1.10D0,1.00D0,1.00D0,1.00D0,1.00D0,
     & 2.20D0,1.80D0,                                           ! 19-36
     &        1.60D0,1.40D0,1.35D0,1.40D0,1.40D0,
     &                        1.40D0,1.35D0,1.35D0,1.35D0,1.35D0,
     &                        1.30D0,1.25D0,1.15D0,1.15D0,1.15D0,1.15D0,
     & 2.35D0,2.00D0,                                           ! 37-54
     &        1.80D0,1.55D0,1.45D0,1.45D0,1.35D0,
     &                        1.30D0,1.35D0,1.40D0,1.60D0,1.55D0,
     &                        1.55D0,1.45D0,1.45D0,1.40D0,1.40D0,1.40D0,
     & 2.6D0,2.15D0,                                           ! 55-86
     &        1.95D0,1.85D0,1.85D0,1.85D0,1.85D0,1.85D0,1.85D0,
     &               1.80D0,1.75D0,1.75D0,1.75D0,1.75D0,1.75D0,1.75D0,
     &        1.75D0,1.55D0,1.45D0,1.35D0,1.35D0,
     &                        1.30D0,1.35D0,1.35D0,1.35D0,1.50D0,
     &                        1.90D0,1.80D0,1.60D0,1.90D0,1.90D0,1.90D0,
     & 2.6D0,2.15D0,                                          ! 87-102
     &        1.95D0,1.80D0,1.75D0,1.75D0,1.75D0,1.75D0,1.75D0,
     &               1.75D0,1.75D0,1.75D0,1.75D0,1.75D0,1.75D0,1.75D0/
*
      Ang2au=1.0D0/angstr
      If (iAtmNr.gt.nAtmNr) Then
        Write(6,*) 'Bragg-Slater: Too high atom number!'
        Write(6,*) 'iAtmNr=',iAtmNr
        Call Quit_OnUserError
      End If
      Bragg_Slater=BS_Radii(iAtmNr)*Ang2au
      Return
      End
