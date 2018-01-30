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
      Subroutine Put_P2MO(P2MO,nP2MO)
      Implicit Real*8 (A-H,O-Z)
#include "SysDef.fh"

      Real*8 P2MO(nP2MO)
      Character*24 Label

      Label='P2MO'
      Call Put_dArray(Label,P2MO,nP2MO)

      Return
      End

      Subroutine Put_P2MOt(P2MO,nP2MO)
      Implicit Real*8 (A-H,O-Z)
#include "SysDef.fh"

      Real*8 P2MO(nP2MO)
      Character*24 Label

      Label='P2MOT'
      Call Put_dArray(Label,P2MO,nP2MO)

      Return
      End
