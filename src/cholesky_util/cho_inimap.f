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
      SUBROUTINE CHO_INIMAP()
      use ChoArr, only: IntMap
C
C     Purpose: initialize integral shell pair calculation map.
C
C     NB!!!!! file is assumed open (restart only)
C
      IMPLICIT NONE
#include "cholesky.fh"

      INTEGER IOPT, NDIM, IADR

      IF (RSTCHO) THEN
         IOPT = 2
         NDIM = SIZE(INTMAP)
         IADR = 0
         CALL IDAFILE(LUMAP,IOPT,INTMAP,NDIM,IADR)
      ELSE
         CALL CHO_IZERO(INTMAP,SIZE(INTMAP))
      END IF

      END
