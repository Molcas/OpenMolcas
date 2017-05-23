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
      SUBROUTINE CHO_IODIAG(DIAG,IOPT)
      IMPLICIT NONE
      REAL*8  DIAG(*)
      INTEGER IOPT

      CHARACTER*7 FNAME
      PARAMETER (FNAME = 'CHODIAG')

      CALL CHO_IODIAG_1(DIAG,IOPT,FNAME)

      END
*
      SUBROUTINE CHO_IOCHODIAG(DIAG,IOPT)
      IMPLICIT NONE
      REAL*8  DIAG(*)
      INTEGER IOPT

      CHARACTER*5 FNAME
      PARAMETER (FNAME = 'CDIAG')

      CALL CHO_IODIAG_1(DIAG,IOPT,FNAME)

      END
*
      SUBROUTINE CHO_IODIAG_1(DIAG,IOPT,FNAME)
C
C     Purpose: write/read a copy of diagonal to disk (1st reduced set).
C              The file is opened and closed here.
C              IOPT=1: write
C              IOPT=2: read
C
#include "implicit.fh"
      DIMENSION     DIAG(*)
      CHARACTER*(*) FNAME
#include "cholesky.fh"

      CHARACTER*12 SECNAM
      PARAMETER (SECNAM = 'CHO_IODIAG_1')

      LUNIT = 7
      CALL DANAME(LUNIT,FNAME)
      IF (IOPT.EQ.1 .OR. IOPT.EQ.2) THEN
         LENGTH = NNBSTRT(1)
         IADR   = 0
         CALL DDAFILE(LUNIT,IOPT,DIAG,LENGTH,IADR)
      ELSE   ! error
         WRITE(LUPRI,*) SECNAM,': IOPT out of bounds: ',IOPT
         CALL CHO_QUIT('Error in '//SECNAM,104)
      END IF
      CALL DACLOS(LUNIT)

      END
