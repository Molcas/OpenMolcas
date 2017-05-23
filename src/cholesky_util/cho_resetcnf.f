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
      SUBROUTINE CHO_RESETCNF
C
C     Purpose: reset configuration of decomposition to that read from
C              restart file. Original configuration will be saved in
C              restart common block.
C
      IMPLICIT NONE
#include "cholesky.fh"

      CALL CHO_DSWAP(THRCOM,XTHRCOM)
      CALL CHO_DSWAP(THRDIAG,XTHRDIAG)
      CALL CHO_DSWAP(DAMP(1),XDAMP(1))
      CALL CHO_DSWAP(DAMP(2),XDAMP(2))
      CALL CHO_DSWAP(SPAN,XSPAN)
      CALL CHO_DSWAP(THRNEG,XTHRNEG)
      CALL CHO_DSWAP(WARNEG,XWARNEG)
      CALL CHO_DSWAP(TOONEG,XTOONEG)
      CALL CHO_LSWAP(SCDIAG,XSCDIAG)

      END
