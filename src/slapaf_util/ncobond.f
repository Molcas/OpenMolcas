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
      Integer Function nCoBond(iAtom,nAtoms,nMax,iTabBonds,
     &                         nBondMax,iTabAtoms)
      Implicit Real*8 (a-h,o-z)
      Integer iTabBonds(3,nBondMax), iTabAtoms(2,0:nMax,nAtoms)
#include "bondtypes.fh"
*                                                                      *
************************************************************************
*                                                                      *
      nCoBond=0
      nn=iTabAtoms(1,0,iAtom)
      Do i = 1, nn
         If(iTabBonds(3,iTabAtoms(2,i,iAtom)).eq.Covalent_Bond)
     &      nCoBond=nCoBond+1
      End Do
      Return
      End
      Integer Function nFgBond(iAtom,nAtoms,nMax,iTabBonds,
     &                         nBondMax,iTabAtoms)
      Implicit Real*8 (a-h,o-z)
      Integer iTabBonds(3,nBondMax), iTabAtoms(2,0:nMax,nAtoms)
#include "bondtypes.fh"
*                                                                      *
************************************************************************
*                                                                      *
      nFgBond=0
      nn=iTabAtoms(1,0,iAtom)
      Do i = 1, nn
         If(iTabBonds(3,iTabAtoms(2,i,iAtom)).eq.Fragments_Bond)
     &      nFgBond=nFgBond+1
      End Do
      Return
      End
