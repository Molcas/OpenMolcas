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
*  Get_Coord_New_All
*
*> @brief
*>   Get the new coordinates from RUNFILE
*> @author R. Lindh
*>
*> @details
*> Place Cartesian coordinates (in a.u.) into array \p Coord_All(3,*).
*>
*> @param[out] Coord_All  Array of coordinates
*> @param[in]  nAtoms_All Number of atoms
************************************************************************
      Subroutine Get_Coord_New_All(Coord_All,nAtoms_All)
      Implicit Real*8 (a-h,o-z)
#include "stdalloc.fh"
      Real*8 Coord_All(3,nAtoms_All)
      Real*8, Dimension (:,:), Allocatable :: CU
      Interface
        Subroutine Get_Coord_New(CU,nAtoms)
        Real*8, Dimension (:,:), Allocatable :: CU
        Integer nAtoms
        End Subroutine
      End Interface
*
      Call Get_nAtoms_All(nAtoms_Allx)
      If (nAtoms_All.ne.nAtoms_Allx) Then
         Write (6,*) 'Get_Coord_All: nAtoms_All.ne.nAtoms_Allx'
         Write (6,*) 'nAtoms_All=',nAtoms_All
         Write (6,*) 'nAtoms_Allx=',nAtoms_Allx
         Call Abend
      End If
      Call Get_Coord_New(CU,nAtoms)
      Call Get_Coord_All_(CU,nAtoms,Coord_All,nAtoms_All)
      Call mma_deallocate(CU)
*
      Return
      End
