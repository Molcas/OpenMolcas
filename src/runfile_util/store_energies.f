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
      Subroutine Store_Energies(nEnergies,Energies,iEnergy)
      Implicit None
      Integer nEnergies,iEnergy
      Real*8 Energies(nEnergies)
*
      Call Put_iScalar('Number of roots',nEnergies)
      Call Put_dArray('Last energies',Energies,nEnergies)
      Call Put_dScalar('Last energy',Energies(iEnergy))
*
      Return
      End
