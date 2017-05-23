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
      SUBROUTINE CONF_ARC_W(IOCC_MIN,IOCC_MAX, NORB,  NEL,IVERTEXW,
     &                        IARCW)
*
* Obtain arcweights for single and double occupied arcs
* from vertex weights
*
* Jeppe Olsen, October 2001
*
#include "implicit.fh"
*. Input
      INTEGER IVERTEXW(NORB+1,NEL+1)
      INTEGER IOCC_MIN(NORB),IOCC_MAX(NORB)
*. Output
      INTEGER IARCW(NORB,NEL,2)
      IZERO = 0
      CALL ISETVC(IARCW,IZERO,2*NORB*NEL)
* IARCW(I,J,K) is weight of arc with occupation K ending at (I,J)
* IARCW(I,J,K) = Sum(J-K < L <= J)   IVERTEXW(I-1,L)
      DO I = 1, NORB
       DO J = 1, NEL
        IF(IOCC_MIN(I).LE.J .AND. J.LE.IOCC_MAX(I)) THEN
          DO K = 1, NEL
            IF(K.EQ.1) IARCW(I,J,K) = IVERTEXW(I-1+1,J+1)
            IF(K.EQ.2.AND.J.GE.2) IARCW(I,J,K) = IVERTEXW(I-1+1,J+1)
     &                                         + IVERTEXW(I-1+1,J-1+1)
          END DO
        END IF
       END DO
      END DO
*
      NTEST = 00
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Arc weights for single occupied arcs '
        CALL IWRTMA(IARCW(1,1,1),NORB,NEL,NORB,NEL)
        WRITE(6,*) ' Arc weights for double occupied arcs '
        CALL IWRTMA(IARCW(1,1,2),NORB,NEL,NORB,NEL)
      END IF
*
      RETURN
      END
