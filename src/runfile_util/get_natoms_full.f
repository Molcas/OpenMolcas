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
*  Get_nAtoms_Full
*
*> @author I. Fdez. Galv&aacute;n
*>
*> @details
*> Get number of all atoms (not only symmetry unique) from RUNFILE,
*> including MM atoms otherwise invisible to gateway/slapaf.
*>
*> @param[out] nAtoms_Full Number of all atoms in the system
************************************************************************
      Subroutine Get_nAtoms_Full(nAtoms_Full)
      Implicit None
      Integer nAtoms_Full,nAtom,nAtMM
      Logical Found
*
      Call Get_nAtoms_All(nAtom)
      Call Qpg_dArray('MMO Coords',Found,nAtMM)
      nAtoms_Full=nAtom+nAtMM/3
*
      Return
      End
