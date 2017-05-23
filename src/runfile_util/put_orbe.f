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
      Subroutine Put_OrbE(OrbE,nOrbE)
      Implicit Real*8 (a-h,o-z)
#include "SysDef.fh"

      Real*8       OrbE(nOrbE)
      Character*24 Label

      Label='OrbE'
      Call Put_dArray(Label,OrbE,nOrbE)

      Return
      End
