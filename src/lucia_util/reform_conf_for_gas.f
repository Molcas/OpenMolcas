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
      SUBROUTINE REFORM_CONF_FOR_GAS(ICONF_GAS,ICONF,IBORB, IBEL,MXPORB,
     &                                  NEL,  IWAY)
C    &                                NORB,NEL,IWAY)
*
* Reform between local and global numbering of
* configuration for given GAS space
*
* IWAY = 1 : Global => Local
* IWAY = 2 : Local => GLobal
*
* Jeppe Olsen, November 2001
*
#include "implicit.fh"
*
      INTEGER ICONF_GAS(MXPORB)
      INTEGER ICONF(*)
*
      IF(IWAY.EQ.1) THEN
        DO IEL = 1, NEL
          ICONF_GAS(IEL) = ICONF(IBEL-1+IEL)  - IBORB + 1
        END DO
      ELSE IF (IWAY.EQ.2) THEN
        DO IEL = 1, NEL
          ICONF(IBEL-1+IEL) = ICONF_GAS(IEL) + IBORB - 1
        END DO
      ELSE
        WRITE(6,*) ' Problem in REFORM_CONF ... , IWAY = ', IWAY
*        STOP       ' Problem in REFORM_CONF ... , IWAY = '
         CALL SYSABENDMSG('lucia_util/reform_conv',
     &                    'Internal error',' ')
      END IF
*
      NTEST = 00
      IF(NTEST.GE.100) THEN
        IF(IWAY.EQ.1) THEN
          WRITE(6,*) ' Global => Local reform of conf '
        ELSE
          WRITE(6,*) ' Local => Global reform of conf '
        END IF
        WRITE(6,*) ' ICONF_GAS : '
        CALL IWRTMA(ICONF_GAS,1,NEL,1,NEL)
        WRITE(6,*) ' Accessed part of ICONF '
        CALL IWRTMA(ICONF,1,IBEL-1+NEL,1,IBEL-1+NEL)
      END IF
*
      RETURN
      END
C       NOP_FOR_CONF(JCONF,NEL)
