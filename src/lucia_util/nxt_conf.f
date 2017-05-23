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
      SUBROUTINE NXT_CONF(ICONF,NEL,NORB,INI,NONEW)
*
* Next configuration of NEL electrons distributed in NORB orbitals
*
* A configuration is stored as the occupied orbitals
* in nonstrict ascending order - two consecutivw orbitals are allowed
* to be identical
* allowing two
*
* IF INI = 1 : Generate initial configuration
*    NONOEW = 1 : No new configuration could be generated
*
* Jeppe Olsen, November 2001
*
#include "implicit.fh"
      INTEGER ICONF(NEL)
*
      NTEST = 00
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Input configuration to NXT_CONF '
        CALL IWRTMA(ICONF,1,NEL,1,NEL)
        WRITE(6,*) ' NEL, NORB = ',  NEL, NORB
      END IF
      IF(INI.EQ.1) THEN
*. Check that NEL electrons can be distributed in NORB orbitals
        IF(NEL.LE.2*NORB) THEN
          NONEW = 0
          N_DOUBLE = NEL/2
          DO I = 1, N_DOUBLE
            ICONF(2*I-1) = I
            ICONF(2*I)   = I
          END DO
          IF(2*N_DOUBLE.NE.NEL) ICONF(2*N_DOUBLE+1) = N_DOUBLE + 1
        ELSE
          NONEW = 1
        END IF
*
      ELSE IF(INI.EQ.0) THEN
*
        IADD  = 1
        IEL = 0
*. Increase orbital number of next electron
 1000   CONTINUE
          IEL = IEL + 1
*. Can orbital number be increased for electron IEL ?
          INCREASE = 0
          IF(IEL.LT.NEL) THEN
            IF(ICONF(IEL).LT.ICONF(IEL+1)-1)  INCREASE = 1
            IF(ICONF(IEL).EQ.ICONF(IEL+1)-1) THEN
*. If ICONF(IEL) is increased, ICONF(IEL) = ICONF(IEL+1), check if this is ok
              IF(IEL.EQ.NEL-1) THEN
                 INCREASE = 1
              ELSE IF (ICONF(IEL+1).NE.ICONF(IEL+2)) THEN
                 INCREASE = 1
              END IF
            END IF
          ELSE
C-jwk new            IF(ICONF(IEL).LT.NORB) THEN
            IF(IEL .EQ. NEL .AND. ICONF(IEL) .LT. NORB) THEN
              INCREASE = 1
            ELSE
*. Nothing more to do
              NONEW = 1
              GOTO 1001
            END IF
          END IF
*
          IF(INCREASE.EQ.1) THEN
*. Increase orbital for elec IEL
            NONEW = 0
            ICONF(IEL) = ICONF(IEL)+1
*. Minimize orbital occupations
            NDOUBLE = (IEL-1)/2
            DO JORB = 1, NDOUBLE
              ICONF(2*JORB-1) = JORB
              ICONF(2*JORB  ) = JORB
            END DO
            IF(2*NDOUBLE.LT.IEL-1) ICONF(IEL-1) = NDOUBLE+1
            IADD = 0
          END IF
        IF(IADD.EQ.1)  GOTO 1000
      END IF
*     ^ End if INI = 0
*
 1001 CONTINUE
*
      IF(NTEST.GE.100) THEN
        IF(NONEW.EQ.1) THEN
          WRITE(6,*) ' No new configurations '
          WRITE(6,*) ' Input configuration '
          CALL IWRTMA(ICONF,1,NEL,1,NEL)
        ELSE
          WRITE(6,*) ' Next configurations '
          CALL IWRTMA(ICONF,1,NEL,1,NEL)
        END IF
      END IF
*
      RETURN
      END
