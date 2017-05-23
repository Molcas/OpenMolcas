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
      SUBROUTINE CONF_GRAPH(IOCC_MIN,IOCC_MAX, NORB,  NEL,IARCW,
     &                        NCONF,   ISCR)
*
* A group of configurations is described by the
* accumulated min and max, IOCC_MIN and IOCC_MAX.
*
* Find arcweights of corresponding graph and total number
* of configurations ( all symmetries)
*
* Jeppe Olsen, Oct. 2001
*
#include "implicit.fh"
*. Input
      INTEGER IOCC_MIN(NORB), IOCC_MAX(NORB)
*. Output
      INTEGER IARCW(NORB,NEL,2)
* IARCW(I,J,K) gives weight of arc ending at vertex (I,J)
* with occupation K (=1,2)
*. Local scratch : Length should be (NORB+1)*(NEL+1)
      INTEGER ISCR(NORB+1,NEL+1)
*. Set up vertex weights
      CALL CONF_VERTEX_W(IOCC_MIN,IOCC_MAX,NORB,NEL,ISCR)
      NCONF = ISCR(NORB+1,NEL+1)
*. Obtain arcweights from vertex weights
C?    WRITE(6,*) ' CONF_GRAPH, NORB, NEL = ', NORB, NEL
      CALL CONF_ARC_W( IOCC_MIN, IOCC_MAX,     NORB,      NEL,     ISCR,
     &                    IARCW)
*
      NTEST = 00
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' IOCMIN and IOCMAX  '
        CALL IWRTMA(IOCC_MIN,1,NORB,1,NORB)
        CALL IWRTMA(IOCC_MAX,1,NORB,1,NORB)
        WRITE(6,*) ' Arcweights for single occupied arcs '
        CALL IWRTMA(IARCW(1,1,1),NORB,NEL,NORB,NEL)
        WRITE(6,*) ' Arcweights for double occupied arcs '
        CALL IWRTMA(IARCW(1,1,2),NORB,NEL,NORB,NEL)
        WRITE(6,*) ' Total number of configurations ', NCONF
      END IF
*
      RETURN
      END
*
