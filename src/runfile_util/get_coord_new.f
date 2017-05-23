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
      Subroutine Get_Coord_New(CN,nAtoms)
************************************************************
*
*   <DOC>
*     <Name>Get\_Coord\_New</Name>
*     <Syntax>Call Get\_Coord\_New(ipCoord,nAtoms)</Syntax>
*     <Arguments>
*       \Argument{ipCoord}{Pointer to Work of the array of the symmetry unique cartesian coordinates of the basis set centers.}{Integer}{out}
*       \Argument{nAtoms}{Number of symmetry unique cartesian coordinates of the basis set centers.}{Integer}{in}
*     </Arguments>
*     <Purpose>To get the updated/new symmetry unique cartesian coordinates of the basis set centers.</Purpose>
*     <Dependencies></Dependencies>
*     <Author>R. Lindh</Author>
*     <Modified_by></Modified_by>
*     <Side_Effects></Side_Effects>
*     <Description>The utility will read the updated/new symmetry unique cartesian coordinates of the basis set centers from the run file.
*     </Description>
*    </DOC>
*
************************************************************
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
