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
      Subroutine Expand_Coor(Coord,nAtoms,W1,iAll_Atom,nSym,iOper)
************************************************************************
*                                                                      *
*     purpose: to generate a list of all atoms                         *
*                                                                      *
************************************************************************
      use Phase_Info
      Implicit Real*8 (A-H,O-Z)
#include "Molcas.fh"
      Real*8 Coord(3,nAtoms)
      Real*8 W1(3,nAtoms*8)
      Integer iGen(3), iCoSet(0:7,0:7), iStab(0:7), iOper(0:nSym-1)
*----------------------------------------------------------------------*
      call dcopy_(nAtoms*3,Coord,1,W1,1)
*----------------------------------------------------------------------*
*     Apply the symmetry operations                                    *
*----------------------------------------------------------------------*
      nGen=0
      If (nSym.eq.2) nGen=1
      If (nSym.eq.4) nGen=2
      If (nSym.eq.8) nGen=3
      If (nGen.ge.1) iGen(1)=iOper(1)
      If (nGen.ge.2) iGen(2)=iOper(2)
      If (nGen.ge.3) iGen(3)=iOper(4)
*
      MaxDCR=0
      iAll_Atom=nAtoms
      Do iAtom = 1, nAtoms
         iChAtom=iChxyz(W1(1,iAtom),iGen,nGen)
         Call Stblz(iChAtom,iOper,nSym,nStab,iStab,MaxDCR,iCoSet)
         nCoSet=nSym/nStab
         XOld=W1(1,iAtom)
         YOld=W1(2,iAtom)
         ZOld=W1(3,iAtom)
*
         Do iCo = 1, nCoSet-1
*
            iAll_Atom = iAll_Atom + 1
            W1(1,iAll_Atom)=XOld*DBLE(iPhase(1,iCoSet(iCo,0)))
            W1(2,iAll_Atom)=YOld*DBLE(iPhase(2,iCoSet(iCo,0)))
            W1(3,iAll_Atom)=ZOld*DBLE(iPhase(3,iCoSet(iCo,0)))
*
         End Do
*
      End Do
*----------------------------------------------------------------------*
*                                                                      *
*----------------------------------------------------------------------*
      Return
      End
