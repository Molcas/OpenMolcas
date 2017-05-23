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
      SUBROUTINE Lucia_Close
      IMPLICIT REAL*8           (A-H,O-Z)
#include "clunit.fh"
*
* Free memory allocated by Lucia
*
      Call FREESTR_GAS()
      Call DeAlloc_Lucia()
*
* Close any files opened by Lucia
*
      CALL DAClos(LUDIA)
      CALL DAClos(LUC)
      CALL DAClos(LUHC)
      CALL DAClos(LUSC1)
      CALL DAClos(LUSC2)
      CALL DAClos(LUSC3)
      CALL DAClos(LUSC34)
      CALL DAClos(LUSC35)
      CALL DAClos(LUSC36)
      CALL DAClos(LUSC37)
      CALL DAClos(LUSC38)
      CALL DAClos(LUSC39)
      CALL DAClos(LUSC40)
      CALL DAClos(LUMOUT)
      RETURN
      END
