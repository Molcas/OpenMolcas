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
      Subroutine Put_Coord_Full(Coord,nAtoms)
      Implicit None
      Integer      nAtoms, nAtoms_All
      Real*8       Coord(3,nAtoms)

      Call Get_nAtoms_All(nAtoms_All)
      Call Put_Coord_New(Coord,nAtoms_All)
      Call Put_dArray('MMO Coords',Coord(1,nAtoms_All+1),
     &                3*(nAtoms-nAtoms_All))

      Return
      End
