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
      SUBROUTINE GETPRINTLEVEL
      IMPLICIT NONE
#include "prgm.fh"
      Integer, External :: iPrintLevel
      Logical, External :: Reduce_Prt

C Global print level
      IPGLOB=iPrintLevel(-1)
C Change print level when in an optimization loop
      IF (Reduce_Prt()) THEN
        IPGLOB=MAX(IPGLOB-USUAL,SILENT)
      END IF

      END
