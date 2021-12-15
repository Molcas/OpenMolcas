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
*  Get_Coord_Full
*
*> @brief
*>   Get coordinates from RUNFILE
*> @author I. Fdez. Galv&aacute;n
*>
*> @details
*> Place Cartesian coordinates (in a.u.) into array \p Coord_Full(3,*).
*> Includes MM atoms otherwise invisible to gateway/slapaf.
*>
*> @param[out] Coord_Full  Array of coordinates
*> @param[in]  nAtoms_Full Number of atoms
************************************************************************
      Subroutine Get_Coord_Full(Coord_Full,nAtoms_Full)
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
        Call Abend
      End If
      Call Get_nAtoms_All(nAtoms_All)
      If (nAtoms_Full.lt.nAtoms_All) Then
        Write (6,*) 'Get_Coord_Full: nAtoms_Full.lt.nAtoms_All'
        Write (6,*) 'nAtoms_Full=',nAtoms_Full
        Write (6,*) 'nAtoms_Fullx=',nAtoms_All
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
