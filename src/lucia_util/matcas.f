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
      SUBROUTINE MATCAS(     CIN,    COUT,   NROWI,   NROWO,  IROWO1,
     &                     NGCOL,    ISCA,  SCASGN)
*
* COUT(IR+IROWO1-1,ISCA(IC)) =
* COUT(IR+IROWO1-1,ISCA(IC)) + CIN(IR,IC)*SCASGN(IC)
* (if IGAT(IC).ne.0)
*
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION CIN(NROWI,*),COUT(NROWO,*)
      INTEGER ISCA(*)
      DIMENSION SCASGN(*)
*
      MAXCOL = 0
      DO 100 IC = 1, NGCOL
        IF(ISCA(IC).NE.0) THEN
          ICEXP = ISCA(IC)
          MAXCOL = MAX(MAXCOL,ICEXP)
          SIGN = SCASGN(IC)
          DO 50 IR = 1,NROWI
            COUT(IR+IROWO1-1,ICEXP) =
     &      COUT(IR+IROWO1-1,ICEXP) + SIGN*CIN(IR,IC)
   50     CONTINUE
        END IF
  100 CONTINUE
*
      NTEST = 0
      IF(NTEST.NE.0) THEN
        WRITE(6,*) ' Output from MATCAS '
        CALL WRTMAT(COUT,NROWO,MAXCOL,NROWO,MAXCOL)
      END IF
*
      RETURN
      END
C                 CALL MATCG(C,CB(ICGOFF),NROW,NIBTC,IBOT,
C                            NKBTC,I1,XI1S)
