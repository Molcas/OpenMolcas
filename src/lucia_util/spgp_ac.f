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
* Copyright (C) 1997, Jeppe Olsen                                      *
************************************************************************
* Some GAS routines
*
* Nomenclature
*
*. Group of strings : set of strings with a given number of orbitals
*                    in a given GASspace
*
*. Supergroup of strings  : product of NGAS groups, i.e. a string with a
*                    given numb er of electrons in each GAS space
*
*. Type of string : Type is defined by the total number of electrons
*                   in the string. A type will therefore in general
*                   consists of several supergroups
*
      SUBROUTINE SPGP_AC( INSPGRP,NINSPGRP,IOUTSPGRP,NOUTSPGRP,  NGAS,
     &                    MXPNGAS,     IAC,ISPGRP_AC,IBASEIN,IBASEOUT)
*
* Annihilation/Creation mapping of strings
*
* Jeppe Olsen, April 1, 1997
*
      IMPLICIT REAL*8(A-H,O-Z)
*. General input : Number of electrons in each gasspace
      INTEGER INSPGRP(MXPNGAS,*),IOUTSPGRP(MXPNGAS,*)
*. Output
      INTEGER ISPGRP_AC(NGAS,*)
*. Check first that supergroups + IAC information is consistent
      NELIN = 0
      NELOUT = 0
      DO IGAS = 1, NGAS
        NELIN  = NELIN + INSPGRP(IGAS,IBASEIN)
        NELOUT = NELOUT + IOUTSPGRP(IGAS,IBASEOUT)
      END DO
      IF(.NOT.((IAC.EQ.1.AND.NELIN.EQ.NELOUT+1).OR.
     &         (IAC.EQ.2.AND.NELIN.EQ.NELOUT-1))) THEN
        WRITE(6,*) ' Inconsistent data provided to SPGP_AC'
        WRITE(6,*) ' NELIN NELOUT IAC=',NELIN,NELOUT,IAC
*        STOP' Inconsistent data provided to SPGRP_AC'
        CALL SYSABENDMSG('lucia_util/spgp_ac','Internal error',' ')
      END IF
*
      DO ISPGRP = IBASEIN, IBASEIN+NINSPGRP-1
        DO IGAS = 1, NGAS
          IF(IAC.EQ.1) THEN
            INSPGRP(IGAS,ISPGRP) = INSPGRP(IGAS,ISPGRP) - 1
          ELSE IF (IAC.EQ.2) THEN
             INSPGRP(IGAS,ISPGRP) = INSPGRP(IGAS,ISPGRP) + 1
          END IF
*. Find corresponding supergroup
          ITO = 0
          DO JSPGRP = IBASEOUT, IBASEOUT+NOUTSPGRP-1
            IAMOKAY = 1
            DO JGAS = 1, NGAS
              IF( INSPGRP(JGAS,ISPGRP).NE.IOUTSPGRP(JGAS,JSPGRP))
     &        IAMOKAY=0
            END DO
            IF(IAMOKAY.EQ.1) ITO = JSPGRP
          END DO
          ISPGRP_AC(IGAS,ISPGRP) = ITO
*. And clean up
          IF(IAC.EQ.1) THEN
            INSPGRP(IGAS,ISPGRP) = INSPGRP(IGAS,ISPGRP) + 1
          ELSE IF (IAC.EQ.2) THEN
             INSPGRP(IGAS,ISPGRP) = INSPGRP(IGAS,ISPGRP) - 1
          END IF
        END DO
      END DO
*
      NTEST = 000
      IF(NTEST.GE.1000) THEN
        WRITE(6,*) ' Input supergroups '
        CALL IWRTMA(INSPGRP(1,IBASEIN),NGAS,NINSPGRP,
     &  MXPNGAS,NINSPGRP)
        WRITE(6,*) ' Output supergroups '
        CALL IWRTMA(IOUTSPGRP(1,IBASEOUT),NGAS,NOUTSPGRP,
     &  MXPNGAS,NOUTSPGRP)
      END IF
*
      IF(NTEST.GE.100) THEN
       WRITE(6,*) ' Output from SPGP_AC '
       WRITE(6,*) ' ===================='
       WRITE(6,*)
       WRITE(6,*) ' IAC = ', IAC
       WRITE(6,*) ' Mapping : '
       CALL IWRTMA(ISPGRP_AC(1,IBASEIN),NGAS,NINSPGRP,NGAS,NINSPGRP)
      END IF
*
      RETURN
      END
