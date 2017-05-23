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
      SUBROUTINE ASIND(IAS,ISYM,ICASE,IP,IQ,IR)
      USE SUPERINDEX
      IMPLICIT REAL*8 (A-H,O-Z)
#include "rasdim.fh"
#include "caspt2.fh"
#include "output.fh"

      GOTO (12,13) ICASE

  12  CONTINUE
      IABABS=IAS+NAGEBES(ISYM)
      IAABS=MAGEB(1,IABABS)
      IBABS=MAGEB(2,IABABS)
      GOTO 1213
  13  CONTINUE
      IABABS=IAS+NAGTBES(ISYM)
      IAABS=MAGTB(1,IABABS)
      IBABS=MAGTB(2,IABABS)
 1213 CONTINUE
      IP=IEXTIS(IAABS)
      IQ=IEXTIS(IBABS)
      IR=0
      RETURN

      END
