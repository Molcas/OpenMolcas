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
*  Put_CMO
*
*> @brief
*>   Write the symmetry blocked MO coefficients on the runfile
*> @author R. Lindh
*>
*> @details
*> The utility will write the symmetry blocked MO coefficients on the runfile.
*>
*> @param[in] CMO  Array of symmetry blocked MO coefficients
*> @param[in] nCMO Number of elements in \p CMO
************************************************************************
      Subroutine Put_CMO(CMO,nCMO)
      Implicit Real*8 (A-H,O-Z)
#include "real.fh"
#include "WrkSpc.fh"

      Real*8    CMO(nCMO)
      Character Label*24

      Label='Last orbitals'
      Call Put_dArray(Label,CMO,nCMO)

      Return
      End
