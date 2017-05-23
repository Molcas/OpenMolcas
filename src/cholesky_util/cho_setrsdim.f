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
      SUBROUTINE CHO_SETRSDIM(NDIMRS,MSYM,MRED,IRED,ILOC)
C
C     Purpose: set reduced set dimension.
C
#include "implicit.fh"
      INTEGER NDIMRS(MSYM,MRED)
#include "cholesky.fh"

      IF (IRED .LE. MAXRED) THEN
         CALL ICOPY(NSYM,NNBSTR(1,ILOC),1,NDIMRS(1,IRED),1)
      END IF

      END
