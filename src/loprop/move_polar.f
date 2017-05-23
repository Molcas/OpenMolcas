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
      Subroutine Move_Polar(Polar,EC,nAtoms,nij,iANr,Bond_Threshold)
*
* Distributes the contributions from the bonds that doesn't fullfill the requiremen
* Bond_Length <= Bond_Threshold*(Bragg_Slater(iAtom)+Bragg_Slater(jAtom)) to the
* two atoms involved in the bond.
*
* Multipole moments are moved in the Move_Prop subroutine!
*
      Implicit Real*8 (a-h,o-z)
#include "real.fh"
      Real*8 EC(3,nij), Polar(6,nij)
      Integer iAnr(nAtoms)
      Logical Bond_OK, Check_Bond
*
      Do iAtom = 2,nAtoms
         ii = iAtom*(iAtom+1)/2
         Do jAtom = 1,iAtom-1
            jj = jAtom*(jAtom+1)/2
            Bond_Ok = Check_Bond(EC(1,ii),EC(1,jj),iANr(iAtom),
     &                           iANr(jAtom),Bond_Threshold)
            ij = iAtom*(iAtom-1)/2+jAtom
            If (.NOT. Bond_Ok) Then
*
* Move half of the bond polarizabilities to each atom
*
               Do iPol = 1,6
                  Polar(iPol,ii) = Polar(iPol,ii) + Half*Polar(iPol,ij)
                  Polar(iPol,jj) = Polar(iPol,jj) + Half*Polar(iPol,ij)
                  Polar(iPol,ij) = Zero
               End Do
            End If
         End Do
      End Do
*
      Return
      End
