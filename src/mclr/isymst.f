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
      FUNCTION ISYMST_MCLR(STRING,NEL)
*
* Master routine for symmetry of string
*
      IMPLICIT REAL*8(A-H,O-Z)
*. General input ( PNTGRP is used )
#include "detdim.fh"
#include "orbinp_mclr.fh"
*. Specific input
      INTEGER STRING(*)
*
      ISYM = 1
      DO 100 IEL = 1, NEL
        JSYM=ISMFTO(STRING(IEL))-1
cVV: this ugly code needed for CRAY SV1
        IVV=ISYM-1
        KVV = IEOR(IVV,JSYM)
        ISYM=KVV+1
100   CONTINUE
      ISYMST_MCLR = ISYM
*
      RETURN
      END
