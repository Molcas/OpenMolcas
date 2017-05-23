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
      Subroutine Compute_T(Z_Tot,T,ZA,RA,nAtoms)
      Implicit Real*8 (a-h,o-z)
#include "real.fh"
      Real*8 T(3), ZA(nAtoms), RA(3,nAtoms)
*                                                                      *
************************************************************************
*                                                                      *
*---- Form the center of nuclear charge
*
      Do iCar = 1, 3
         Tmp = Zero
         Do iAtom = 1, nAtoms
            Tmp = Tmp + ZA(iAtom)*RA(iCar,iAtom)
         End Do
         T(iCar) = Tmp/Z_Tot
      End Do
*                                                                      *
************************************************************************
*                                                                      *
      Return
      End
