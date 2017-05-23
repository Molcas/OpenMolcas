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
      SUBROUTINE MAX_OPEN_ORB(MAXOP,IOCLS,NGAS,NOCLS,NOBPT)
*
* Max number of open orbitals in occupation classes
*
* Jeppe Olsen, November 2001
*
#include "implicit.fh"
#include "csfbas.fh"
*. Input
      INTEGER IOCLS(NGAS,NOCLS)
      INTEGER NOBPT(NGAS)
*
      MAXOP = 0
C?    WRITE(6,*) ' NOCLS, NGAS = ', NOCLS, NGAS
      DO JOCLS = 1, NOCLS
        MAXOP_J = 0
        DO IGAS = 1, NGAS
          NEL = IOCLS(IGAS,JOCLS)
          NORB = NOBPT(IGAS)
C?        WRITE(6,*) ' IGAS, NEL, NORB = ', IGAS, NEL, NORB
          MAXOP_IGAS = MIN(NEL,2*NORB-NEL)
          MAXOP_J = MAXOP_J + MAXOP_IGAS
        END DO
        MAXOP = MAX(MAXOP,MAXOP_J)
      END DO
      maxop_lucia=maxop
*
      NTEST = 00
      IF(NTEST.GE.100) THEN
        WRITE(6,*)
     &  ' Max number of unpaired orbitals = ', MAXOP
      END IF
*
      RETURN
      END
