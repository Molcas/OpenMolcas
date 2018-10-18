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
* Copyright (C) 2018, Ignacio Fdez. Galvan                             *
************************************************************************
*  Get_Mass
*
*> @brief
*>   Get (symmetry-unique) atomic masses from RUNFILE
*> @author Ignacio Fdez. Galn&aacute;n
*>
*> @details
*> Place atomic masses (in a.u.) into array \p Mass_All(*).
*>
*> @param[out] Mass_All   Array of masses
*> @param[in]  nAtoms_All Number of atoms
************************************************************************
      Subroutine Get_Mass(Mass,nAtoms)
      Implicit None
#include "stdalloc.fh"
      Real*8 Mass(nAtoms)
      Integer nAtoms,mAtoms,nCent,i
      Integer, Dimension (:), Allocatable :: AtoB
      Real*8, Dimension (:), Allocatable :: CentMass
      Logical Found

      Call Get_iScalar('Unique atoms',mAtoms)
      If (mAtoms.ne.nAtoms) Then
        Write (6,*) 'Get_Mass: mAtoms.ne.nAtoms'
        Write (6,*) 'mAtoms=',mAtoms
        Write (6,*) 'nAtoms=',nAtoms
        Call QTrace()
        Call Abend()
      End If
      Call mma_allocate(AtoB,nAtoms)
      Call Get_iArray('Atom -> Basis',AtoB,nAtoms)
      Call Qpg_dArray('Isotopes',Found,nCent)
      If (.not.Found) Then
        Write (6,*) 'Get_Mass: Isotopes array not found'
        Call QTrace()
        Call Abend()
      End If
      Call mma_allocate(CentMass,nCent)
      Call Get_dArray('Isotopes',CentMass,nCent)
      Do i=1,nAtoms
        Mass(i) = CentMass(AtoB(i))
      End Do
      Call mma_deallocate(CentMass)
      Call mma_deallocate(AtoB)

      Return
      End Subroutine Get_Mass
