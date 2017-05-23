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
      Subroutine Move_Prop(rMP,EC,lMax,nElem,nAtoms,nPert,
     &                     nij,iANr,Bond_Threshold)
*
* Distributes the contributions from the bonds that doesn't fulfill the requirement
* Bond_Length <= Bond_Threshold*(Bragg_Slater(iAtom)+Bragg_Slater(jAtom) to the
* two atoms involved in the bond.
*
* Polarizabilities are moved in the Move_Polar subroutine!
*
      Implicit Real*8 (A-H,O-Z)
#include "real.fh"
      Real*8 rMP(nij,0:nElem-1,0:nPert-1), EC(3,nij)
      Integer iAnr(nAtoms)
      Logical Bond_OK, Check_Bond
*
      Do iAtom = 2,nAtoms
         ii = iAtom*(iAtom+1)/2
         Do jAtom = 1,iAtom-1
            jj = jAtom*(jAtom+1)/2
            Bond_Ok = Check_Bond(EC(1,ii),EC(1,jj),iANr(iAtom),
     &                           iANr(jAtom),Bond_Threshold)
            If (.NOT. Bond_OK) Then
               ij = iAtom*(iAtom-1)/2+jAtom
               Do iPert = 0,nPert-1
*
* First move half of the bond properties to iAtom
*
                  Do iElem = 0,nElem-1
                     rMP(ij,iElem,iPert) = rMP(ij,iElem,iPert) * Half
                  End Do
*
                  Call ReExpand(rMP(1,0,iPert),nij,nElem,EC(1,ij),
     &                          EC(1,ii),ij,lMax)
*
                  Do iElem = 0,nElem-1
                     rMP(ii,iElem,iPert) = rMP(ii,iElem,iPert)
     &                                   + rMP(ij,iElem,iPert)
                  End Do
*
* Then move the other half of the bond properties to jAtom
*
                  Call ReExpand(rMP(1,0,iPert),nij,nElem,EC(1,ii),
     &                          EC(1,jj),ij,lMax)
*
                  Do iElem = 0,nElem-1
                     rMP(jj,iElem,iPert) = rMP(jj,iElem,iPert)
     &                                   + rMP(ij,iElem,iPert)
                  End Do
*
* Set local properties to zero
*
                  Call dCopy_(nElem,Zero,0,rMP(ij,0,iPert),nij)
               End Do
            End If
         End Do
      End Do
*
      Return
      End
