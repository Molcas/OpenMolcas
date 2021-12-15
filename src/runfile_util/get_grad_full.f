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
*  Get_Grad_Full
*
*> @brief
*>   Get gradient from RUNFILE
*> @author I. Fdez. Galv&aacute;n
*>
*> @details
*> Place Cartesian gradient (in a.u.) into array \p Grad_Full(3,*).
*> Includes MM atoms otherwise invisible to gateway/slapaf.
*>
*> @param[out] Grad_Full  Array of gradient
*> @param[in]  nAtoms_All Number of atoms
************************************************************************
      Subroutine Get_Grad_Full(Grad_Full,nAtoms_Full)
      Implicit None
      Integer nAtoms_Full, nAtoms_Fullx, nAtoms_All, nGrad, nGradMM
      Real*8 Grad_Full(3,nAtoms_Full)
      Logical Found
*
      Call Get_nAtoms_Full(nAtoms_Fullx)
      If (nAtoms_Full.ne.nAtoms_Fullx) Then
        Write (6,*) 'Get_Grad_Full: nAtoms_Full.ne.nAtoms_Fullx'
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
      Call Qpg_dArray('GRAD',Found,nGrad)
      If(.not.Found .or. nGrad.eq.0) Then
        Write (6,*) 'Get_Grad_Full: Did not find GRAD'
        Call Abend
      End If
      Call Get_dArray('GRAD',Grad_Full,nGrad)
      Call Qpg_dArray('MMO Grad',Found,nGradMM)
      If (Found) Then
        Call Get_dArray('MMO Grad',Grad_Full(1,nAtoms_All+1),nGradMM)
      End If
*
      Return
      End
