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
      SUBROUTINE INFO_CONF_LIST(NCONF_PER_OPEN,
     &                           MAXOP,NEL,LENGTH_LIST,NCONF_TOT,IB_REO,
     &                           IB_OCC)
*
* Info on  configuration list form NCONF_PER_OPEN
*
* Jeppe Olsen, November 2001
*
#include "implicit.fh"
*. Input
      INTEGER NCONF_PER_OPEN(MAXOP+1)
*. Output :
*    Offset for configuration  with given number of
*    open orbitals in of configuration reordering
      INTEGER IB_REO(MAXOP+1)
*    Offset for configuration  with given number of
*    open orbitals in list of configuration occupations
      INTEGER IB_OCC(MAXOP+1)
*
      LENGTH = 0
      JB_REO = 1
      JB_OCC = 1
C?    WRITE(6,*) ' MAXOP, NEL = ', MAXOP, NEL
      DO NOPEN = 0, MAXOP
        IB_REO(NOPEN+1) =  JB_REO
        IB_OCC(NOPEN+1) =  JB_OCC
        IF(MOD(NEL-NOPEN,2).EQ.0) THEN
          NOCOB = NOPEN + (NEL - NOPEN)/2
          JB_OCC = JB_OCC + NOCOB*NCONF_PER_OPEN(NOPEN+1)
          JB_REO = JB_REO + NCONF_PER_OPEN(NOPEN+1)
        END IF
      END DO
*
      LENGTH_LIST = JB_OCC-1
      NCONF_TOT = JB_REO - 1
*
      NTEST = 00
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' NCONF_PER_OPEN list '
        CALL IWRTMA(NCONF_PER_OPEN,1,MAXOP+1,1,MAXOP+1)
        WRITE(6,*) ' Length of configuration list :', LENGTH_LIST
        WRITE(6,*) ' Total number of configurations : ', NCONF_TOT
      END IF
*
      RETURN
      END
