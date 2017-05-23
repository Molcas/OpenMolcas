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
      SUBROUTINE MATCG(     CIN,    COUT,   NROWI,   NROWO,  IROWI1,
     &                    NGCOL,    IGAT,  GATSGN)
*
* Gather columns of CIN with phase
*
* COUT(IR,IC) = GATSGN(IC)*CIN(IR+IROWI1-1,IGAT(IC)) if IGAT(IC) .ne.0
* COUT(IR,IC) = 0                           if IGAT(IC) .ne.0
*
      IMPLICIT REAL*8(A-H,O-Z)
      INTEGER IGAT(*)
      DIMENSION GATSGN(*)
      DIMENSION CIN(NROWI,*),COUT(NROWO,*)
*
C?    write(6,*) ' MATCG NROWI,NROWO,IROWI1,NGCOL '
C?    write(6,*)         NROWI,NROWO,IROWI1,NGCOL
      DO 100 IG = 1, NGCOL
C?      write(6,*) ' igat,sign ',IGAT(IG),GATSGN(IG)
        IF(IGAT(IG).EQ.0) THEN
          DO 20 IR = 1, NROWO
            COUT(IR,IG)=0.0D0
   20     CONTINUE
        ELSE
         IGFRM = IGAT(IG)
         SIGN = GATSGN(IG)
         DO 30 IR = 1, NROWO
           COUT(IR,IG) = SIGN*CIN(IROWI1-1+IR,IGFRM)
   30    CONTINUE
        END IF
  100 CONTINUE
*
      NTEST = 0
      IF(NTEST.NE.0) THEN
        WRITE(6,*) ' Column gathered matrix '
        CALL WRTMAT(COUT,NROWO,NGCOL,NROWO,NGCOL)
      END IF
*
      RETURN
      END
