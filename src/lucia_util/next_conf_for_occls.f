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
      SUBROUTINE NEXT_CONF_FOR_OCCLS            (ICONF,
     &                                           IOCCLS,
     &                                           NGAS,
     &                                           NOBPT,
     &                                           INI,
*
     &                                           NONEW)
*
* Obtain next configuration for occupation class
*
* Jeppe Olsen, Nov. 2001
*
#include "implicit.fh"
#include "mxpdim.fh"
*
*. Input
*
*. Number of electrons per gas space
      INTEGER IOCCLS(NGAS)
*. Number of orbitals per gasspace
      INTEGER NOBPT(NGAS)
*. Input and output
      INTEGER ICONF(*)
*. Local scratch
      INTEGER IBORB(MXPNGAS), ICONF_GAS(MXPORB)
      INTEGER IBEL(MXPNGAS)
*
      NTEST = 000
*. Total number of electrons
      NEL = IELSUM(IOCCLS,NGAS)
C?    WRITE(6,*) ' NEXT_CONF ... NEL, NGAS = ', NEL, NGAS
*. Offset for orbitals and electrons
      DO IGAS = 1, NGAS
        IF(IGAS.EQ.1) THEN
          IBORB(IGAS) = 1
          IBEL(IGAS)  = 1
        ELSE
          IBORB(IGAS) = IBORB(IGAS-1)+NOBPT(IGAS-1)
          IBEL(IGAS)  = IBEL(IGAS-1) + IOCCLS(IGAS-1)
        END IF
      END DO
*
      NONEW = 1

      IF(INI.EQ.1) THEN
*
*. Initial configuration
*
        NONEW = 0
        INI_L = 1
        NONEW_L = 0
        DO IGAS = 1, NGAS
*. Initial configuration for this GASSPACE
          NEL_GAS = IOCCLS(IGAS)
          NORB_GAS = NOBPT(IGAS)
C?        WRITE(6,*) ' IGAS, NEL_GAS, NORB_GAS = ',
C?   &                 IGAS, NEL_GAS, NORB_GAS
          CALL NXT_CONF(ICONF_GAS,NEL_GAS,NORB_GAS,INI_L,NONEW_L)
          IF(NONEW_L.EQ.1) THEN
             NONEW = 1
             GOTO 1001
          ELSE
             JBEL   = IBEL(IGAS)
             JBORB  = IBORB(IGAS)
             JEL = IOCCLS(IGAS)
             JORB  = NOBPT(IGAS)
             CALL REFORM_CONF_FOR_GAS             (ICONF_GAS,
     &                                             ICONF,
     &                                             JBORB,
     &                                             JBEL,MXPORB,JEL,2)
C    &            (ICONF_GAS,ICONF,JBORB,JBEL,JORB,JEL,2)
C                 (ICONF_GAS,ICONF,IBORB,IBEL,NORB,NEL,IWAY)
          END IF
        END DO
*
        IF(NTEST.GE.1000) THEN
          WRITE(6,*) ' Initial configuration '
          CALL IWRTMA(ICONF,1,NEL,1,NEL)
        END IF
*
      ELSE
*
*. Next configuration
*
*. Loop over GAS spaces and find first GASspace where a new configuration
*. could be obtained
        DO IGAS = 1, NGAS
C?        WRITE(6,*) ' IGAS = ', IGAS
*. Remove the offsets for this space
          JBEL   = IBEL(IGAS)
          JBORB  = IBORB(IGAS)
          JEL = IOCCLS(IGAS)
          JORB = NOBPT(IGAS)
C?        WRITE(6,*) ' JBEL, JBORB, JEL, JORB = ',
C?   &                 JBEL, JBORB, JEL, JORB
          CALL REFORM_CONF_FOR_GAS          (ICONF_GAS,
     &                                       ICONF,
     &                                       JBORB,JBEL,MXPORB,JEL, 1)
C    &         (ICONF_GAS,ICONF,JBORB,JBEL,JORB,JEL,1)
*. Generate next configuration for this space
          INI_L = 0
          CALL NXT_CONF(ICONF_GAS,JEL,JORB,INI_L,NONEW_L)
          IF(NONEW_L.EQ.0) THEN
            NONEW = 0
*. Configuration in space IGAS, was increased. Copy this and reset configurations
*. in previous gasspaces to initial form
            CALL REFORM_CONF_FOR_GAS            (ICONF_GAS,
     &                                           ICONF,
     &                                           JBORB,
     &                                           JBEL,
     &                                           MXPORB,
*
     &                                           JEL,
     &                                             2)
C    &           (ICONF_GAS,ICONF,JBORB,JBEL,JORB,JEL,2)
*
            DO JGAS = 1, IGAS-1
              JBEL   = IBEL(JGAS)
              JBORB  = IBORB(JGAS)
              JEL = IOCCLS(JGAS)
              JORB = NOBPT(JGAS)
              CALL REFORM_CONF_FOR_GAS              (ICONF_GAS,
     &                                               ICONF,
     &                                               JBORB,
     &                                               JBEL,MXPORB,JEL,1)
C    &             (ICONF_GAS,ICONF,JBORB,JBEL,JORB,JEL,1)
              INI_L = 1
              CALL NXT_CONF(ICONF_GAS,JEL,JORB,INI_L,NONEW_L)
C                                     NEL_GAS,NORB_GAS,INI_L,NONEW_L)
              CALL REFORM_CONF_FOR_GAS              (ICONF_GAS,
     &                                               ICONF,
     &                                               JBORB,
     &                                               JBEL,MXPORB,JEL,2)
C    &             (ICONF_GAS,ICONF,JBORB,JBEL,JORB,JEL,2)
            END DO
*. Get out of the loop
            GOTO 1001
          END IF
        END DO
*       ^ End of loop over gasspaces
      END IF
*     ^ End if swith between initialization/next
 1001 CONTINUE
*
      IF(NTEST.GE.100) THEN
        IF(NONEW.EQ.1) THEN
          WRITE(6,*) ' No new configuration '
        ELSE
          WRITE(6,*) ' New configuration '
          CALL IWRTMA(ICONF,1,NEL,1,NEL)
        END IF
      END IF
*
      RETURN
      END
