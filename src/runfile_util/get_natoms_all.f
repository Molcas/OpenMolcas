************************************************************************
* This file is part of OpenMolcas.                                     *
*                                                                      *
* OpenMolcas is free software; you can redistribute it and/or modify   *
* it under the terms of the GNU Lesser General Public License, v. 2.1. *
* OpenMolcas is distributed in the hope that it will be useful, but it *
* is provided "as is" and without any express or implied warranties.   *
* For more details see the full text of the license in the file        *
* LICENSE or in <http://www.gnu.org/licenses/>.                        *
*                                                                      *
* Copyright (C) Roland Lindh                                           *
************************************************************************
*  Get_nAtoms_All
*
*> @author R. Lindh
*>
*> @details
*> Get number of all atoms (not only symmetry unique) from RUNFILE.
*>
*> @param[out] nAtoms_All Number of all atoms in the molecule
************************************************************************
      Subroutine Get_nAtoms_All(nAtoms_All)
      Implicit Real*8 (a-h,o-z)
#include "real.fh"
#include "WrkSpc.fh"
*
      Call Get_iScalar('Unique atoms',nAtoms)
      Call Allocate_Work(ipCoord,3*nAtoms)
      Call Get_dArray('Unique Coordinates',Work(ipCoord),3*nAtoms)
      Call Get_nAtoms_All_(Work(ipCoord),nAtoms,nAtoms_all)
      Call Free_Work(ipCoord)
*
*     Write (*,*) 'nAtoms_All=',nAtoms_All
*     Write (*,*) 'Exit Get_nAtoms_All'
      Return
      End
      Subroutine Get_nAtoms_All_(Coord_Unique_Atoms,nUnique_Atoms,
     &                           nAll_Atoms)
************************************************************
      use Symmetry_Info, only: nIrrep, iOper, Symmetry_Info_Get
      Implicit Real*8 (a-h,o-z)
#include "real.fh"
      Integer iGen(3), iCoSet(0:7)
      Real*8  Coord_Unique_Atoms(3,nUnique_Atoms)
      Integer, Save:: Active=0
*     Write (*,*) 'Enter Get_nAtoms_All_'
*                                                                      *
************************************************************************
*                                                                      *
      If (Active.eq.0) Then
         Call Symmetry_Info_Get()
         Active=1
      End If
*     Write (*,*) 'Get_nAtoms_All_: nIrrep',nIrrep
*     Write (*,*) 'Get_nAtoms_All_: iOper',(iOper(i),i=0,nIrrep-1)
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
*     Write (*,*) 'nGen=',nGen
*     Write (*,*) 'iGen=',iGen
*                                                                      *
************************************************************************
*                                                                      *
*     Compute total number of centers.
*
      iAll_Atom=0
      Do iUnique_Atom = 1, nUnique_Atoms
*        Write (*,*) 'iUnique_Atom=',iUnique_Atom
*
         iChAtom=iChxyz(Coord_Unique_Atoms(1,iUnique_Atom),iGen,nGen)
*        Write (*,*) 'iChAtom=',iChAtom
         Call CoSet(iCoSet,nCoSet,iChAtom,iOper,nIrrep)
*        Write (*,*) 'nCoSet=',nCoSet
*        Write (*,*) 'iCoSet=',(iCoSet(i),i=0,nCoset-1)
*
         iAll_Atom = iAll_Atom + nCoSet
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
