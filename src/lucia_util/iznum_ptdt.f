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
* Copyright (C) 2001, Jeppe Olsen                                      *
************************************************************************
      FUNCTION IZNUM_PTDT(IAB,NOPEN,NALPHA,Z,NEWORD,IREORD)
*
* Adress of prototype determinant IAB
* alpha occupation is used to define lex address
*
* Jeppe Olsen, Dec. 2001
*
#include "implicit.fh"
      INTEGER Z(NOPEN,NALPHA)
      DIMENSION IAB(*),NEWORD(*)
*
      NTEST = 00
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' IZNUM_PTDT, NOPEN, NALPHA = ', NOPEN,NALPHA
        WRITE(6,*) ' Input Z- matrix '
        CALL IWRTMA(Z,NOPEN,NALPHA,NOPEN,NALPHA)
      END IF
*
      IZ = 1
      IALPHA = 0
      DO I = 1,NOPEN
        IF(IAB(I).GT.0) THEN
          IALPHA = IALPHA + 1
          IZ = IZ + Z(I,IALPHA)
        END IF
      END DO
*
C?    WRITE(6,*) ' IZ = ', IZ
      IF(IREORD.EQ.0) THEN
        IZNUM_PTDT = IZ
      ELSE
        IZNUM_PTDT = NEWORD(IZ)
      END IF
*
      IF ( NTEST .GE. 100 ) THEN
        WRITE(6,*) ' Output from IZNUM_PTDT '
        WRITE(6,*) ' Prototype determinant '
        CALL IWRTMA(IAB,1,NOPEN,1,NOPEN)
        WRITE(6,*) ' Address = ', IZNUM_PTDT
      END IF
*
      RETURN
      END
