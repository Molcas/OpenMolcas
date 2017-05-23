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
      SUBROUTINE CHO_UNINI()
C
C     Purpose: zero unit numbers.
C
      IMPLICIT NONE
#include "cholesky.fh"

      LURED  = 0
      CALL CHO_IZERO(LUCHO,NSYM)
      LURST  = 0
      LUMAP  = 0

      END
