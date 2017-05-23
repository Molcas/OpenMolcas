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
      SUBROUTINE CONF_VERTEX_W(IOCC_MIN,IOCC_MAX,NORB,NEL,IVERTEXW)
*
* Obtain vertex weights for configuration graph
*
* Jeppe Olsen, October 2001
*
#include "implicit.fh"
*. Input
      INTEGER IOCC_MIN(NORB),IOCC_MAX(NORB)
*. Output
      INTEGER IVERTEXW(NORB+1,NEL+1)
*
      IZERO = 0
C?    WRITE(6,*) ' CONF_VERTEX : NORB, NEL = ', NORB, NEL
      CALL ISETVC(IVERTEXW,IZERO,(NORB+1)*(NEL+1))
*
      IVERTEXW(0+1,0+1) = 1
      DO IORB = 1, NORB
C$$$ Jesper      DO IEL  = 0, NEL
        DO IEL  = IOCC_MIN(IORB), IOCC_MAX(IORB)
C$$$ Jesper        IF(IOCC_MIN(IORB).LE.IEL.AND.IEL.LE.IOCC_MAX(IORB)) THEN
*
          IF(IEL.EQ. 0 )
     &    IVERTEXW(IORB+1,IEL+1) = IVERTEXW(IORB-1+1,IEL+1)
*
          IF(IEL.EQ. 1 )
     &    IVERTEXW(IORB+1,IEL+1) = IVERTEXW(IORB-1+1,IEL+1)
     &                           + IVERTEXW(IORB-1+1,IEL+1-1)
*
          IF(IEL.GE. 2 )
     &    IVERTEXW(IORB+1,IEL+1) = IVERTEXW(IORB-1+1,IEL+1)
     &                           + IVERTEXW(IORB-1+1,IEL+1-1)
     &                           + IVERTEXW(IORB-1+1,IEL+1-2)
*
C$$$ Jesper        END IF
        END DO
      END DO
* Check whether a configuration has enough singly occupied orbitals
C$$$ Jesper ???  CALL REDUCE_VERTEX_WEIGHTS(IVERTEXW(NORB+1,NEL+1), IVERTEXW,
C$$$ Jesper ??? &                           NEL, NORB, IOCC_MIN, IOCC_MAX)
*
      NTEST = 00
      IF(NTEST.GE.100) THEN
        WRITE(6,*)  ' Vertex weights as an (NORB+1)*(NEL+1) matrix '
        CALL IWRTMA(IVERTEXW,NORB+1,NEL+1,NORB+1,NEL+1)
      END IF
*
      RETURN
      END
*
*
