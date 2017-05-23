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
      SUBROUTINE MKNSM
C     PUPROSE: CREATE THE SYMMETRY INDEX VECTOR
C
      IMPLICIT REAL*8 (A-H,O-Z)
C
#include "rasdim.fh"
#include "rasscf.fh"
#include "general.fh"
#include "gas.fh"
C
      NLEV=0
      DO IGAS=1,NGAS
        DO ISYM=1,NSYM
          NSTA=NLEV+1
          NLEV=NLEV+NGSSH(IGAS,ISYM)
          DO LEV=NSTA,NLEV
            NSM(LEV)=ISYM
          END DO
        END DO
      END DO
      RETURN
      END
