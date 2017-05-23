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
      SUBROUTINE CHO_OPFVEC(ISYM,IOPT)
C
C     Purpose: open/close files for full storage vectors, sym. ISYM.
C
#include "implicit.fh"
#include "cholesky.fh"
#include "choreo.fh"

      CHARACTER*10 SECNAM
      PARAMETER (SECNAM = 'CHO_OPFVEC')

      CHARACTER*6 FNAME

      MULD2H(I,J)=IEOR(I-1,J-1)+1

      IF (IOPT .EQ. 0) THEN
         DO ISYMA = 1,NSYM
            DO ISYMB = 1,ISYMA
               LUFV(ISYMA,ISYMB) = -1
               LUFV(ISYMB,ISYMA) = -1
            END DO
         END DO
      ELSE IF (IOPT .EQ. 1) THEN
         DO ISYMB = 1,NSYM
            ISYMA = MULD2H(ISYMB,ISYM)
            IF (ISYMA .GE. ISYMB) THEN
               WRITE(FNAME,'(A4,I1,I1)') REONAM,ISYMA,ISYMB
               LUNIT = 7
               CALL DANAME_MF_WA(LUNIT,FNAME)
               LUFV(ISYMA,ISYMB) = LUNIT
               LUFV(ISYMB,ISYMA) = LUNIT
            END IF
         END DO
      ELSE IF (IOPT .EQ. 2) THEN
         DO ISYMB = 1,NSYM
            ISYMA = MULD2H(ISYMB,ISYM)
            IF (ISYMA .GE. ISYMB) THEN
               LUNIT = LUFV(ISYMA,ISYMB)
               CALL DACLOS(LUNIT)
               LUFV(ISYMA,ISYMB) = -1
               LUFV(ISYMB,ISYMA) = -1
            END IF
         END DO
      ELSE
         CALL CHO_QUIT('IOPT error in '//SECNAM,104)
      END IF

      END
