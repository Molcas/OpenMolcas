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
      Subroutine Move_Xhole(rX,EC,nAtoms,nij,iANr,Bond_Threshold)

*
* Distribute the XHole-integrals if the bond lengths does not fulfill:
* Bond_Length <= Bond_Threshold*(Bragg_Slater(iAtom)+Bragg_Slater(jAtom) to the
* two atoms involved in the bond.
* Ripped from move_prop.
*

      Implicit Real*8 (A-H,O-Z)
#include "real.fh"
      Dimension rX(nij), EC(3,nij)
      Dimension iAnr(nAtoms)
      Logical Bond_OK, Check_Bond

      Do iAtom = 2,nAtoms
         ii = iAtom*(iAtom+1)/2
         Do jAtom = 1,iAtom-1
            jj = jAtom*(jAtom+1)/2
            Bond_Ok = Check_Bond(EC(1,ii),EC(1,jj),iANr(iAtom),
     &                           iANr(jAtom),Bond_Threshold)
            If (.NOT. Bond_OK) Then
               ij = iAtom*(iAtom-1)/2+jAtom
*
* First move half of the bond properties to iAtom
*
               rX(ij) = rX(ij) * Half
               rX(ii) = rX(ii) + rX(ij)
*
* Then move the other half of the bond properties to jAtom
*
               rX(jj) = rX(jj) + rX(ij)
*
* Set local properties to zero
*
               rX(ij) = Zero
            End If
         End Do
      End Do

      Return
      End
