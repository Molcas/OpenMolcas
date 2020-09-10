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
*  Get_Coord_All
*
*> @brief
*>   Get coordinates from RUNFILE
*> @author R. Lindh
*>
*> @details
*> Place Cartesian coordinates (in a.u.) into array \p Coord_All(3,*).
*>
*> @param[out] Coord_All  Array of coordinates
*> @param[in]  nAtoms_All Number of atoms
************************************************************************
      Subroutine Get_Coord_All(Coord_All,nAtoms_All)
      Implicit Real*8 (a-h,o-z)
#include "real.fh"
#include "stdalloc.fh"
      Real*8 Coord_All(3,nAtoms_All)
      Real*8, Dimension (:,:), Allocatable :: CU

      Call Get_nAtoms_All(nAtoms_Allx)
      If (nAtoms_All.ne.nAtoms_Allx) Then
         Write (6,*) 'Get_Coord_All: nAtoms_All.ne.nAtoms_Allx'
         Write (6,*) 'nAtoms_All=',nAtoms_All
         Write (6,*) 'nAtoms_Allx=',nAtoms_Allx
         Call QTrace
         Call Abend
      End If
      Call Get_iScalar('Unique atoms',nAtoms)
      Call mma_allocate(CU,3,nAtoms)
      Call Get_dArray('Unique Coordinates',CU,3*nAtoms)
      Call Get_Coord_All_(CU,nAtoms,Coord_All,nAtoms_All)
      Call mma_deallocate(CU)
*
      Return
      End
      Subroutine Get_Coord_All_(Coord_Unique,nUnique_Atoms,
     &                          Coord_All,nAll_Atoms)
      use Symmetry_Info, only: nIrrep, iOper, Symmetry_Info_Get
      Implicit Real*8 (a-h,o-z)
#include "real.fh"
      Integer iGen(3), iCoSet(0:7,0:7), iStab(0:7)
      Real*8  Coord_Unique(3,nUnique_Atoms),
     &        Coord_All(3,nAll_Atoms)
      integer, Save:: Active=0
*     Write (*,*) 'Enter Get_Coord_All_'
*                                                                      *
************************************************************************
*                                                                      *
      if(Active.eq.0) then
       Call Symmetry_Info_Get()
       Active=1
      endif
*     Write (*,*) 'Get_Coord_All_: nIrrep=',nIrrep
*     Write (*,*) 'Get_Coord_All_: iOper=',(iOper(i),i=0,nIrrep-1)
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
*     Generate list of all coordinates, index arrays, etc.
*
      iAll_Atom=0
      MaxDCR=0
      Do iUnique_Atom = 1, nUnique_Atoms
*        Write (*,*) 'iUnique_Atom=',iUnique_Atom
*
         iChAtom=iChxyz(Coord_Unique(1,iUnique_Atom),iGen,nGen)
*        Write (*,*) 'iChAtom=',iChAtom
         Call Stblz(iChAtom,nStab,iStab,MaxDCR,iCoSet)
         nCoSet=nIrrep/nStab
*        Write (*,*) 'In Get_Coord_All'
*        Write (*,*) 'nCoset=',nCoset
*        Write (*,*) 'iCoset=',(iCoset(i,0),i=0,nCoset-1)
*
         Do iCo = 0, nCoSet-1
*           Write (*,*) 'In Get_Coord_All'
*           Write (*,*) 'iCo,iCoSet(iCo,0)=',iCo,iCoSet(iCo,0)
            iAll_Atom = iAll_Atom + 1
            Call OA(iCoSet(iCo,0),Coord_Unique(1:3,iUnique_Atom),
     &              Coord_All(1:3,iAll_Atom))
         End Do
*
      End Do
*
*     Call RecPrt('Coord_Unique',' ',Coord_Unique,3,nUnique_Atoms)
*     Call RecPrt('Coord_All',' ',Coord_All,3,nAll_Atoms)
*                                                                      *
************************************************************************
*                                                                      *
      Return
      End
