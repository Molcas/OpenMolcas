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
* Copyright (C) 2015, Ignacio Fdez. Galvan                             *
************************************************************************
      Subroutine Get_Coord_Full(Coord_Full,nAtoms_Full)
************************************************************
*
*   <DOC>
*     <Name>Get\_Coord\_Full</Name>
*     <Syntax>Call Get\_Coord\_Full(Coord\_Full,nAtoms\_Full)</Syntax>
*     <Arguments>
*       \Argument{Coord\_Full}{Array of coordinates}{Real*8 (3,nAtoms\_All)}{out}
*       \Argument{nAtoms\_Full}{Number of atoms}{Integer}{in}
*     </Arguments>
*     <Purpose>Get Coordinates from RUNFILE</Purpose>
*     <Dependencies>Get\_Coord\_All, Get\_nAtoms\_Full</Dependencies>
*     <Author>I. Fdez. Galvan</Author>
*     <Modified_by></Modified_by>
*     <Side_Effects></Side_Effects>
*     <Description>
*       Place cartesian coordinates (in a.u. into array Coord\_Full(3,*).
*       Includes MM atoms otherwise invisible to gateway/slapaf.
*     </Description>
*    </DOC>
*
************************************************************
      Implicit None
      Integer nAtoms_Full, nAtoms_Fullx, nAtoms_All, nCoordMM
      Real*8 Coord_Full(3,nAtoms_Full)
      Logical Found
*
      Call Get_nAtoms_Full(nAtoms_Fullx)
      If (nAtoms_Full.ne.nAtoms_Fullx) Then
        Write (6,*) 'Get_Coord_Full: nAtoms_Full.ne.nAtoms_Fullx'
        Write (6,*) 'nAtoms_Full=',nAtoms_Full
        Write (6,*) 'nAtoms_Fullx=',nAtoms_Fullx
        Call QTrace
        Call Abend
      End If
      Call Get_nAtoms_All(nAtoms_All)
      If (nAtoms_Full.lt.nAtoms_All) Then
        Write (6,*) 'Get_Coord_Full: nAtoms_Full.lt.nAtoms_All'
        Write (6,*) 'nAtoms_Full=',nAtoms_Full
        Write (6,*) 'nAtoms_Fullx=',nAtoms_All
        Call QTrace
        Call Abend
      End If
      Call Get_Coord_All(Coord_Full,nAtoms_All)
      Call Qpg_dArray('MMO Coords',Found,nCoordMM)
      If (Found) Then
         Call Get_dArray('MMO Coords',Coord_Full(1,nAtoms_All+1),
     &                   nCoordMM)
      End If
*
      Return
      End
