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
      SUBROUTINE EXCIND(IAS,INS,ISYM,ICASE,IP,IQ,IR,IS)
      IMPLICIT REAL*8 (A-H,O-Z)

#include "rasdim.fh"
#include "caspt2.fh"
#include "output.fh"
#include "WrkSpc.fh"
      EXTERNAL ASIND

CPAM99 New call sequence for ASIND
      CALL ASIND(IAS,ISYM,ICASE,IP1,IQ1,IR1)
      CALL NSIND(INS,ISYM,ICASE,IP2,IQ2,IR2)

      IP=IQ1
      IQ=IQ2
      IR=IP1
      IS=IP2
      RETURN
      END
