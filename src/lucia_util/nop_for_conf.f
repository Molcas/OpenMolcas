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
      FUNCTION NOP_FOR_CONF(ICONF,NEL)
*
* A configuration is given as a nonstrict ascending sequence of occupied
* occupied orbitals. Find number of double occupied orbitals
*
* Jeppe Olsen, Nov. 2001
*
#include "implicit.fh"
      INTEGER ICONF(NEL)
*. Loop over electrons
      NOPEN = 0
      IEL = 1
 1000 CONTINUE
        IF(IEL.LT.NEL) THEN
          IF(ICONF(IEL).NE.ICONF(IEL+1)) THEN
           NOPEN = NOPEN + 1
           IEL = IEL + 1
          ELSE IF (ICONF(IEL).EQ.ICONF(IEL+1) ) THEN
           IEL = IEL + 2
          END IF
        END IF
*
        IF(IEL.EQ.NEL) THEN
*. The last orbital is not identical to any later orbitals so
         NOPEN = NOPEN+1
         IEL = IEL + 1
        END IF
      IF(IEL.LT.NEL) GOTO 1000
*
      NOP_FOR_CONF = NOPEN
*
      NTEST = 00
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Configuration '
        CALL IWRTMA(ICONF,1,NEL,1,NEL)
        WRITE(6,*) ' Number of open orbitals = ', NOP_FOR_CONF
      END IF
*
      RETURN
      END
