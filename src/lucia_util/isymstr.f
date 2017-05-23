************************************************************************
* This file is part of OpenMolcas.                                     *
*                                                                      *
* OpenMolcas is free software; you can redistribute it and/or modify   *
* it under the terms of the GNU Lesser General Public License, v. 2.1. *
* OpenMolcas is distributed in the hope that it will be useful, but it *
* is provided "as is" and without any express or implied warranties.   *
* For more details see the full text of the license in the file        *
* LICENSE or in <http://www.gnu.org/licenses/>.                        *
*                                                                      *
* Copyright (C) 1998, Jeppe Olsen                                      *
************************************************************************
      FUNCTION ISYMSTR(ISYM,NSTR)
*
* Symmetry of product of NSTR string symmetries
*
* works currently only for D2H and subgroups
*
* Jeppe Olsen, 1998
*
      IMPLICIT REAL*8(A-H,O-Z)
*. General input
#include "multd2h.fh"
*. Input
      INTEGER ISYM(*)
*
      IF(NSTR.EQ.0) THEN
        IISYM = 1
      ELSE
        IISYM = ISYM(1)
        DO JSTR = 2, NSTR
           IISYM = MULTD2H(IISYM,ISYM(JSTR))
        END DO
      END IF
*
      ISYMSTR = IISYM
*
      RETURN
      END
