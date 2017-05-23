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
      INTEGER FUNCTION CHO_IRANGE(I,BIN,NBIN,LOWER)
C
C     Purpose: return range to which I belongs.
C              It is assumed that BIN is sorted
C              in ascending order so that:
C
C                I > BIN(NBIN)   ==> CHO_IRANGE = NBIN
C              Else
C                I > BIN(NBIN-1) ==> CHO_IRANGE = NBIN-1
C              Else
C                I > BIN(NBIN-2) ==> CHO_IRANGE = NBIN-2
C                ..
C                ..
C                ..
C              Else
C                I < BIN(1)      ==> CHO_IRANGE = 0
C
C              If NBIN <= 0      ==> CHO_IRANGE = -1
C
C              In case of degeneracies, the smallest bin is
C              returned if flag LOWER is set, i.e.,
C
C              I > BIN(IBIN) = BIN(IBIN-1)
C                            = ...
C                            = BIN(IBIN-N) ==> CHO_IRANGE = IBIN-N
C
C              else the largest bin (IBIN in the case above).
C
      IMPLICIT NONE
      INTEGER I, NBIN, BIN(*)
      LOGICAL LOWER
      INTEGER IBIN, JBIN, JBIN1

      CHO_IRANGE = 0

      IF (NBIN .LT. 1) THEN
         CHO_IRANGE = -1
      ELSE
         IF (LOWER) THEN
            DO IBIN = NBIN,1,-1
               IF (I .GT. BIN(IBIN)) THEN
                  CHO_IRANGE = IBIN
                  JBIN1 = IBIN - 1
                  DO JBIN = JBIN1,1,-1
                     IF (BIN(JBIN) .EQ. BIN(IBIN)) THEN
                        CHO_IRANGE = JBIN
                     ELSE
                        RETURN
                     END IF
                  END DO
                  RETURN
               END IF
            END DO
         ELSE
            CHO_IRANGE = 1
            DO IBIN = NBIN,2,-1
               IF (I .GT. BIN(IBIN)) THEN
                  CHO_IRANGE = IBIN
                  RETURN
               END IF
            END DO
         END IF
      END IF

      END
