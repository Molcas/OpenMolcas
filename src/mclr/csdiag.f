************************************************************************
* This file is part of OpenMolcas.                                     *
*                                                                      *
* OpenMolcas is free software; you can redistribute it and/or modify   *
* it under the terms of the GNU Lesser General Public License, v. 2.1. *
* OpenMolcas is distributed in the hope that it will be useful, but it *
* is provided "as is" and without any express or implied warranties.   *
* For more details see the full text of the license in the file        *
* LICENSE or in <http://www.gnu.org/licenses/>.                        *
************************************************************************
      SUBROUTINE CSDIAG(CSFDIA,DETDIA,NCNFTP,NTYP,
     &                   ICTSDT,NDTFTP,NCSFTP,IFLAG ,
     &                   NCNFCN,ICNFOK,IPRNT)
*
*.. obtain average CI diagonal elements and store in
*   CSFDIA as CSF diagonal
*
*
*
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION CSFDIA(*),DETDIA(*)
      DIMENSION NCNFTP(NTYP),NDTFTP(NTYP),NCSFTP(NTYP)
      DIMENSION ICTSDT(*),ICNFOK(*)
*
      ICSOFF = 1
      IDTOFF = 1
      JCNABS = 0
      DO ITYP = 1, NTYP
        IDET = NDTFTP(ITYP)
        ICSF = NCSFTP(ITYP)
        ICNF = NCNFTP(ITYP)
        DO  JCNF = 1, ICNF
          JCNABS = JCNABS + 1
          EAVER = 0.0D0
          DO  JDET = 1, IDET
            EAVER = EAVER +DETDIA(ABS(ICTSDT(IDTOFF-1+JDET) ))
          End Do
          IF( IDET .NE. 0 )EAVER = EAVER/DBLE(IDET)
          CALL SETVEC(CSFDIA(ICSOFF),EAVER,ICSF)
          ICSOFF = ICSOFF + ICSF
          IDTOFF = IDTOFF + IDET
       End Do
      End Do
      RETURN
c Avoid unused argument warnings
      IF (.FALSE.) THEN
        CALL Unused_integer(IFLAG)
        CALL Unused_integer(NCNFCN)
        CALL Unused_integer_array(ICNFOK)
        CALL Unused_integer(IPRNT)
      END IF
      END
