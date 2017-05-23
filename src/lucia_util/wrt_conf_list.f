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
      SUBROUTINE WRT_CONF_LIST(ICONF,NCONF_FOR_OPEN,MAXOP,NCONF,NELEC)
*
* Write list of configurations, given in packed form
*
* Jeppe Olsen, November 2001
*
#include "implicit.fh"
*
      INTEGER ICONF(*), NCONF_FOR_OPEN(MAXOP+1)
*
      IB = 1
      DO IOPEN = 0, MAXOP
        NCONF_OP = NCONF_FOR_OPEN(IOPEN+1)
        IF(NCONF_OP.NE.0) THEN
          WRITE(6,*) ' Number of configurations with ', IOPEN,
     &               ' open orbitals is ', NCONF_OP
*
          NOCC_ORB = IOPEN + (NELEC-IOPEN)/2
          DO JCONF = 1, NCONF_OP
            CALL IWRTMA(ICONF(IB),1,NOCC_ORB,1,NOCC_ORB)
            IB = IB + NOCC_ORB
          END DO
        END IF
      END DO
*
      RETURN
c Avoid unused argument warnings
      IF (.FALSE.) CALL Unused_integer(NCONF)
      END
