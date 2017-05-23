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
      Subroutine Put_Coord_New(Coord,nAtoms)
************************************************************
*
*   <DOC>
*     <Name>Put\_Coord\_New</Name>
*     <Syntax>Call Put\_Coord\_New(Coord,nAtoms)</Syntax>
*     <Arguments>
*       \Argument{Coord}{Array of the updated/new symmetry unique cartesian coordinates of the basis set centers.}
*                       {Real*8 (3,nAtoms)}{in}
*       \Argument{nAtoms}{Number of symmetry unique basis set centers.}{Integer}{in}
*     </Arguments>
*     <Purpose>To write the  updated/new symmetry unique cartesian coordinates of the basis set centers on the run file.</Purpose>
*     <Dependencies></Dependencies>
*     <Author>R. Lindh</Author>
*     <Modified_by></Modified_by>
*     <Side_Effects></Side_Effects>
*     <Description>The utility will write  updated/new symmetry unique cartesian coordinates of the basis set centers on the run file.
*     </Description>
*    </DOC>
*
************************************************************
      Implicit Real*8 (a-h,o-z)
#include "real.fh"

      Real*8       Coord(3,nAtoms)
      Character*24 Label

      Label='GeoNew'
      Call Put_dArray(Label,Coord,3*nAtoms)

      Return
      End
