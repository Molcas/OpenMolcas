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
* Copyright (C) 1995, Jeppe Olsen                                      *
************************************************************************
      SUBROUTINE GTSPGP(IEL,ISPGP,IWAY)
*
*
* Relation between number of electrons in AS1, AS2 ... and
* supergoup number
*
* IWAY = 1 :
* Get ISPGP : Supergroup of strings
*             with IEL(*)  electrons in each AS
* IWAY = 2 :
* GET IEL(*)  : Number of electrons in each AS for supergroup ISPGP
*
*
* Jeppe Olsen, Another lonely night in Lund
*               GAS version July 1995
*
      IMPLICIT REAL*8 (A-H,O-Z)
*. Generel input
#include "mxpdim.fh"
#include "WrkSpc.fh"
#include "strbas.fh"
#include "stinf.fh"
#include "cgas.fh"
#include "gasstr.fh"
*. input(IWAY = 2 ), output (IWAY = 1 )
      INTEGER IEL(*)
*
      IF(IWAY.EQ.1) THEN
*. Occupation => Number
        ISPGP = -1
        DO JSPGP = 1, NTSPGP
          IF(ISPGP.EQ.-1) THEN
            IEQUAL = 1
            DO IGAS = 1, NGAS
              IF(NELFSPGP(IGAS,JSPGP).NE.IEL(IGAS))  IEQUAL= 0
            END DO
            IF(IEQUAL.EQ.1) ISPGP = JSPGP
          END IF
        END DO
      ELSE IF (IWAY .EQ. 2 ) THEN
*. Number => Occupation
        DO IGAS = 1, NGAS
         IEL(IGAS) = NELFSPGP(IGAS,ISPGP)
        END DO
      END IF
*
      NTEST  = 000
      IF(NTEST .GE. 100 ) THEN
        WRITE(6,*) ' Output from GTSPGP '
        WRITE(6,*)
     &   ' IWAY ISPGP IEL ', IWAY,ISPGP,(IEL(IGAS),IGAS = 1, NGAS)
      END IF
*
      RETURN
      END
