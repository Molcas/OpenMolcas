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
* Copyright (C) 2003, Jeppe Olsen                                      *
*               2003, Jesper Wisborg Krogh                             *
************************************************************************
      SUBROUTINE MXMNOC_OCCLS(  MINEL,  MAXEL, NORBTP,NORBFTP, NELFTP,
     &                          MINOP, NTESTG)
*
* Construct accumulated MAX and MIN arrays for an occupation class
*
* MINOP ( Smallest allowed number of open orbitals) added
* April2, 2003, JO (modified by JWK, April - June 2003)
*
      IMPLICIT REAL*8           ( A-H,O-Z)
*. Output
      DIMENSION  MINEL(*),MAXEL(*)
*. Input
      INTEGER NORBFTP(*),NELFTP(*)
*. Local scratch added April 2, 2003
*
#include "mxpdim.fh"
C     INTEGER MINOP_ORB(MXPORB), MINOP_GAS(MXPNGAS), MAXOP_GAS(MXPNGAS)
      INTEGER MINOP_GAS(MXPNGAS), MAXOP_GAS(MXPNGAS)
*
      NTESTL = 00
      NTEST = MAX(NTESTG,NTESTL)
*
      IF(NTEST.GE.100) THEN
        WRITE(6,*)
        WRITE(6,*) ' ============'
        WRITE(6,*) ' MXMNOC_OCCLS'
        WRITE(6,*) ' ============'
        WRITE(6,*)
        WRITE(6,*) ' MINOP  = ', MINOP
        WRITE(6,*) ' NORBTP = ', NORBTP
        WRITE(6,*) ' NORBFTP : '
        CALL IWRTMA(NORBFTP,1,NORBTP,1,NORBTP)
      END IF
*. Well
      NGAS = NORBTP
*
*. Largest number of unpaired electrons in each gas space
*
      DO IGAS = 1, NGAS
        MAXOP_GAS(IGAS) = MIN(NELFTP(IGAS),2*NORBFTP(IGAS)-NELFTP(IGAS))
      END DO
*
*. Smallest number of electrons in each GAS space
*
*. 1 : Just based on number of electrons in each space
      DO IGAS = 1, NGAS
        IF(MOD(NELFTP(IGAS),2).EQ.1) THEN
          MINOP_GAS(IGAS) = 1
        ELSE
          MINOP_GAS(IGAS) = 0
        END IF
      END DO
*. 2 : the total number of open orbitals should be MINOP, this puts
*. also a constraint on the number of open orbitals
*
*. The largest number of open orbitals, all spaces
      MAXOP_T = IELSUM(MAXOP_GAS,NGAS)
      DO IGAS = 1, NGAS
*. Max number of open orbitals in all spaces except IGAS
       MAXOP_EXL = MAXOP_T - MAXOP_GAS(IGAS)
       MINOP_GAS(IGAS) = MAX(MINOP_GAS(IGAS),MINOP-MAXOP_EXL)
       IF (MOD(NELFTP(IGAS)-MINOP_GAS(IGAS),2) .EQ. 1) THEN
          MINOP_GAS(IGAS) = MINOP_GAS(IGAS) + 1
       ENDIF
      END DO
*. We now have the min and max number of open shells per occls,
*. Find the corresponding min and max number accumulated electrons,
*
* The Max occupation is obtained by occupying in max in the
* first orbitals
* The Min occupation is obtained by occopying max in the
* last orbitals.
*
      NEL_INI = 0
      IBORB   = 1
      DO IGAS = 1, NGAS
        NELEC = NELFTP(IGAS)
* PAM2009: This looks like a bug. MAX_DOUBLE can go negative.
*        MAX_DOUBLE = (NELEC-MINOP_GAS(IGAS))/2
* Replace with:
        MAX_DOUBLE = MAX(0,(NELEC-MINOP_GAS(IGAS))/2)
*
* If you are in a situation with no electrons to spare
*
        IF (NELEC .EQ. 0) THEN
           DO IORB = 1,NORBFTP(IGAS)
              IF (IORB+IBORB-1 .EQ. 1) THEN
                 MINEL(IORB+IBORB-1) = 0
                 MAXEL(IORB+IBORB-1) = 0
              ELSE
                 MINEL(IORB+IBORB-1) = MINEL(IORB+IBORB-2)
                 MAXEL(IORB+IBORB-1) = MAXEL(IORB+IBORB-2)
              END IF
           END DO
           GOTO 10
        END IF
*
* The min number of electrons
*
*. Doubly occupy the last MAX_DOUBLE orbitals
C Start Jesper !!!
        IF (NORBFTP(IGAS)-MAX_DOUBLE .LE. 0
     &        .AND. MINOP_GAS(IGAS) .GT. 0) CALL Abend
C End Jesper !!!
        IORB_START = MAX(1,NORBFTP(IGAS)-MAX_DOUBLE)
        DO IORB = IORB_START,NORBFTP(IGAS)
           MINEL(IORB+IBORB-1) =
     &           NEL_INI + NELEC - 2*(NORBFTP(IGAS)-IORB)
C?        write(6,*) ' 1 IORB+IBORB-1, MINEL() ',
C?   &    IORB+IBORB-1,  MINEL(IORB+IBORB-1)
        END DO
*. Singly occupy
        DO IORB = NORBFTP(IGAS)-MAX_DOUBLE-1,1,-1
           MINEL(IORB+IBORB-1) = MAX(NEL_INI,MINEL(IORB+IBORB-1+1)-1)
C?        write(6,*) ' 2 IORB+IBORB-1, MINEL() ',
C?   &    IORB+IBORB-1,  MINEL(IORB+IBORB-1)
        END DO
*
*. The max number of electrons
*
       DO IORB = 1, MAX_DOUBLE
         MAXEL(IORB+IBORB-1) = NEL_INI + 2*IORB
       END DO
       DO IORB = MAX_DOUBLE+1, NORBFTP(IGAS)
         IF (IORB+IBORB-1 .EQ. 1) THEN
           MAXEL(IORB+IBORB-1) = 1
         ELSE
           MAXEL(IORB+IBORB-1)=MIN(NEL_INI+NELEC,MAXEL(IORB+IBORB-2)+1)
         ENDIF
       END DO
  10   CONTINUE
       NEL_INI = NEL_INI + NELFTP(IGAS)
       IBORB = IBORB + NORBFTP(IGAS)
      END DO
*
      IF( NTEST .GE. 100 ) THEN
        NORB = IELSUM(NORBFTP,NORBTP)
        WRITE(6,*) ' MINEL : '
        CALL IWRTMA(MINEL,1,NORB,1,NORB)
        WRITE(6,*) ' MAXEL : '
        CALL IWRTMA(MAXEL,1,NORB,1,NORB)
      END IF
*
      RETURN
      END
* routines for spinadapting CC amplitudes
*
