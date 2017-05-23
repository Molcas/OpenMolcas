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
      SUBROUTINE EXTRT_MS_OPEN_OB(IDET_OC,IDET_MS,IDET_OPEN_MS,NEL)
*
* A determinant IDET_OC, IDET_MS is given. Extract spinprojections
* for open shells
*
* Jeppe Olsen, December 2001
*
#include "implicit.fh"
*.input
      INTEGER IDET_OC(NEL),IDET_MS(NEL)
*. Output
      INTEGER IDET_OPEN_MS(*)
*
      IEL = 1
      IOPEN = 0
*. Loop over electrons
 1000 CONTINUE
       IF(IEL.LT.NEL) THEN
         IF(IDET_OC(IEL).NE.IDET_OC(IEL+1)) THEN
*. Single occupied orbital so
            IOPEN = IOPEN + 1
            IDET_OPEN_MS(IOPEN) = IDET_MS(IEL)
            IEL = IEL + 1
          ELSE
            IEL = IEL + 2
          END IF
       ELSE
*. Last electron was not identical to previous, so
*. neccessarily single occupied.
          IOPEN = IOPEN + 1
          IDET_OPEN_MS(IOPEN) = IDET_MS(IEL)
          IEL = IEL + 1
       END IF
      IF(IEL.LE.NEL) GOTO 1000
*
      NTEST = 00
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Input det, occ and ms '
        CALL IWRTMA(IDET_OC,1,NEL,1,NEL)
        CALL IWRTMA(IDET_MS,1,NEL,1,NEL)
        WRITE(6,*) ' Number of open orbitals = ', IOPEN
        WRITE(6,*) ' Output det : ms of open orbitals '
        CALL IWRTMA(IDET_OPEN_MS,1,IOPEN,1,IOPEN)
      END IF
*
      RETURN
      END
