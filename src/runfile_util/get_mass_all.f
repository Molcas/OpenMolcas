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
*               2018, Ignacio Fdez. Galvan                             *
************************************************************************
*  Get_Mass_All
*
*> @brief
*>   Get atomic masses from RUNFILE
*> @author Ignacio Fdez. Galn&aacute;n
*>
*> @details
*> Place atomic masses (in a.u.) into array \p Mass_All(*).
*>
*> @param[out] Mass_All   Array of masses
*> @param[in]  nAtoms_All Number of atoms
************************************************************************
      Subroutine Get_Mass_All(Mass_All,nAtoms_All)
      Implicit None
#include "stdalloc.fh"
      Real*8 Mass_All(nAtoms_All)
      Integer nAtoms_All, nAtoms_Allx, nAtoms
      Real*8, Dimension (:), Allocatable :: Mass
      Real*8, Dimension (:,:), Allocatable :: CU
      Integer i,j,nSym,nGen,MaxDCR,iOper(0:7),iChAtom,iCo,nCoSet,nStab
      Integer iGen(3),iCoSet(0:7,0:7),iStab(0:7)
      Integer, External :: iChxyz

*     Obtain symmetry-unique masses
      Call Get_nAtoms_All(nAtoms_Allx)
      If (nAtoms_All.ne.nAtoms_Allx) Then
         Write (6,*) 'Get_Coord_All: nAtoms_All.ne.nAtoms_Allx'
         Write (6,*) 'nAtoms_All=',nAtoms_All
         Write (6,*) 'nAtoms_Allx=',nAtoms_Allx
         Call QTrace
         Call Abend
      End If
      Call Get_iScalar('Unique atoms',nAtoms)
      Call mma_allocate(Mass,nAtoms)
      Call Get_Mass(Mass,nAtoms)

*     Replicate masses
      Call mma_allocate(CU,3,nAtoms)
      Call Get_dArray('Unique Coordinates',CU,3*nAtoms)
      Call Get_iScalar('nSym',nSym)
      Call Get_iArray('Symmetry operations',iOper,nSym)
      nGen=0
      If (nSym.eq.2) nGen=1
      If (nSym.eq.4) nGen=2
      If (nSym.eq.8) nGen=3
      If (nGen.ge.1) iGen(1)=iOper(1)
      If (nGen.ge.2) iGen(2)=iOper(2)
      If (nGen.eq.3) iGen(3)=iOper(4)
      MaxDCR=0
      j=0
      Do i=1,nAtoms
        iChAtom=iChxyz(CU(1,i),iGen,nGen)
        Call Stblz(iChAtom,iOper,nSym,nStab,iStab,MaxDCR,iCoSet)
        nCoSet=nSym/nStab
        Do iCo=0,nCoSet-1
          j=j+1
          Mass_all(j)=Mass(i)
        End Do
      End Do
      Call mma_deallocate(CU)
      Call mma_deallocate(Mass)

      Return
      End Subroutine Get_Mass_All
