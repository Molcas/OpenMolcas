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
      SUBROUTINE TRPAD3(MAT,FACTOR,NDIM)
*
*  MAT(I,J) = MAT(I,J) + FACTOR*MAT(J,I)
*
*. With some considerations of effective cache use for large
*  matrices
*
      IMPLICIT REAL*8           (A-H,O-Z)
      REAL*8            MAT(NDIM,NDIM)
              FAC2 = 1.0D0 - FACTOR**2
*
C     IWAY = 1
      IWAY = 2
      IF(IWAY.EQ.1) THEN
*
*. No blocking
*
*. Lower half
        DO  J = 1, NDIM
          DO  I = J, NDIM
            MAT(I,J) =MAT(I,J) + FACTOR * MAT(J,I)
          END DO
        END DO
*. Upper half
        IF( ABS(FACTOR) .NE. 1.0D0 ) THEN
          FAC2 = 1.0D0 - FACTOR**2
          DO I = 1, NDIM
            DO J = 1, I - 1
              MAT(J,I) = FACTOR*MAT(I,J ) + FAC2 * MAT(J,I)
            END DO
          END DO
        ELSE IF(FACTOR .EQ. 1.0D0) THEN
          DO I = 1, NDIM
            DO J = 1, I - 1
              MAT(J,I) = MAT(I,J )
            END DO
          END DO
        ELSE IF(FACTOR .EQ. -1.0D0) THEN
          DO I = 1, NDIM
            DO J = 1, I - 1
              MAT(J,I) =-MAT(I,J )
            END DO
          END DO
        END IF
      ELSE IF(IWAY .EQ. 2 ) THEN
*. Simple blocking of matrix
        LBLK = 40
        NBLK = NDIM/LBLK
        IF(NBLK*LBLK.LT.NDIM) NBLK = NBLK + 1
        IOFF = 1-LBLK
C?      write(6,*) 'NBLK ',nblk
        DO IBLK = 1, NBLK
          IF(IBLK.EQ.-1) write(6,*) 'IBLK = ',IBLK
          IOFF = IOFF + LBLK
          IEND = MIN(IOFF+LBLK-1,NDIM)
          JOFF = 1 - LBLK
          DO JBLK = 1, IBLK
            JOFF = JOFF + LBLK
            JEND = MIN(JOFF+LBLK-1,NDIM)
*. Lower half
            DO  I = IOFF,IEND
              IF(IBLK.EQ.JBLK) JEND = I
              DO J = JOFF,JEND
                MAT(I,J) = MAT(I,J) + FACTOR*MAT(J,I)
              END DO
            END DO
*. Upper half
            IF( ABS(FACTOR) .NE. 1.0D0 ) THEN
              FAC2 = 1.0D0 - FACTOR**2
              DO I = IOFF, IEND
                IF(IBLK.EQ.JBLK) JEND = I
                DO J = JOFF, JEND
                  MAT(J,I) = FACTOR*MAT(I,J ) + FAC2 * MAT(J,I)
                 END DO
               END DO
            ELSE IF(FACTOR .EQ. 1.0D0) THEN
              DO I = IOFF, IEND
                IF(IBLK.EQ.JBLK) JEND = I -1
                DO J = JOFF, JEND
                  MAT(J,I) = MAT(I,J )
                END DO
              END DO
            ELSE IF(FACTOR .EQ. -1.0D0) THEN
              DO I = IOFF, IEND
                IF(IBLK.EQ.JBLK) JEND = I
                DO J = JOFF, JEND
                  MAT(J,I) = - MAT(I,J )
                END DO
              END DO
            END IF
*. ENd of loop over blocks
          END DO
        END DO
*. End of IWAY branching
      END IF
      RETURN
      END
