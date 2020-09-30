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
* Copyright (C) 2013, Thomas Bondo Pedersen                            *
************************************************************************
      Subroutine RPA_Setup()
C
C     Thomas Bondo Pedersen (CTCC,UiO), July 2013.
C
C     Set up RPA calculation:
C        - process input
C        - read reference orbitals
C
      Implicit None
#include "rpa_config.fh"

      Character*9 SecNam
      Parameter (SecNam='RPA_Setup')

      ! Register entry

      ! Define data in common blocks (dummy values).
      Call RPA_SetInc()

      ! Get and check integral representation from Runfile.
      Call RPA_SetIntegralRepresentation()
      Call RPA_CheckIntegralRepresentation()

      ! Pick up data from Runfile.
      Call RPA_RdRun()

      ! Process input.
      Call RPA_RdInp()

      ! Read orbitals and orbital energies.
      Call RPA_RdOrb()

      ! Postprocessing and print
      Call RPA_PPInp()

      ! Add info for testing
      Call RPA_Setup_Add_Info()

      ! Register exit

      End
