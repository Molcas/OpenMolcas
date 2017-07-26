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
*  Get_Coord_New
*
*> @brief
*>   Get the updated/new symmetry unique Cartesian coordinates of the basis set centers
*> @author R. Lindh
*>
*> @details
*> The utility will read the updated/new symmetry unique Cartesian coordinates of the basis set centers from the runfile.
*>
*> @param[out] ipCoord Pointer to \c Work of the array of the symmetry unique Cartesian coordinates of the basis set centers
*> @param[in]  nAtoms  Number of symmetry unique Cartesian coordinates of the basis set centers
************************************************************************
      Subroutine Get_Coord_New(CN,nAtoms)
      Implicit Real*8 (a-h,o-z)
      Real*8, Dimension(:,:), Allocatable :: CN
#include "stdalloc.fh"

      Character*24 Label
      Logical      Found

      Label='GeoNew'
      Call qpg_dArray(Label,Found,nAtoms3)
      nAtoms=nAtoms3/3
      If(.not.Found .or. nAtoms3.eq.0) Return
      Call mma_allocate(CN,3,nAtoms)
      Call Get_dArray(Label,CN,nAtoms3)

      Return
      End
