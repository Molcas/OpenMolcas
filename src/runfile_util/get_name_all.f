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
      Subroutine Get_Name_All(Element)
      Implicit Real*8 (a-h,o-z)
#include "real.fh"
#include "WrkSpc.fh"
#include "Molcas.fh"
      Character*2 Element(*), Element_Uniqe(MxAtom)

      Call Get_iScalar('Unique atoms',nAtoms)
      Call Allocate_Work(ipCoord,3*nAtoms)
      Call Get_dArray('Unique Coordinates',Work(ipCoord),3*nAtoms)
      Call Get_Name(Element_Uniqe)
      Call Get_Name_All_(Work(ipCoord),nAtoms,nAtoms_all,
     &                   Element_Uniqe,Element)
      Call Free_Work(ipCoord)
*
      Return
      End
      Subroutine Get_Name_All_(Coord_Unique_Atoms,nUnique_Atoms,
     &                           nAll_Atoms,Element_Uniqe,Element)
      use Symmetry_Info, only: nIrrep, iOper, Symmetry_Info_Get
      Implicit Real*8 (a-h,o-z)
#include "real.fh"
      Integer iGen(3), iCoSet(0:7)
      Real*8  Coord_Unique_Atoms(3,nUnique_Atoms)
      Character*(*) Element(*), Element_Uniqe(nUnique_Atoms)
      Integer, Save:: Active=0
*                                                                      *
************************************************************************
*                                                                      *
      If (Active.eq.0) Then
         Call Symmetry_Info_Get()
         Active=1
      End If
*                                                                      *
************************************************************************
*                                                                      *
      nGen=0
      If (nIrrep.eq.2) nGen=1
      If (nIrrep.eq.4) nGen=2
      If (nIrrep.eq.8) nGen=3
      If (nGen.ge.1) iGen(1)=iOper(1)
      If (nGen.ge.2) iGen(2)=iOper(2)
      If (nGen.eq.3) iGen(3)=iOper(4)
*                                                                      *
************************************************************************
*                                                                      *
*     Compute total number of centers.
*
      iAll_Atom=0
      Do iUnique_Atom = 1, nUnique_Atoms
*
         iChAtom=iChxyz(Coord_Unique_Atoms(1,iUnique_Atom),iGen,nGen)
         Call CoSet(iCoSet,nCoSet,iChAtom,iOper,nIrrep)
*
         Do i = 1, nCoSet
            iAll_Atom = iAll_Atom + 1
            Element(iAll_Atom)=Element_Uniqe(iUnique_Atom)
         End Do
*
      End Do
*
      nAll_Atoms=iAll_Atom
*                                                                      *
************************************************************************
*                                                                      *
*     Write (*,*) 'Exit Get_nAtoms_All_'
      Return
      End
