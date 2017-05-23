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
      ! Copy-transpose matrix A to B using blocks to optimize cache
      ! usage or extended BLAS functionality if present.
      ! B := transpose(A)

      ! double precision version
      SUBROUTINE DTRANS(NROWS,NCOLS,A,LDA,B,LDB)
      IMPLICIT NONE
      ! arguments
      INTEGER, INTENT(IN) :: NROWS, NCOLS
      INTEGER, INTENT(IN) :: LDA, LDB
      REAL*8, INTENT(IN)  :: A(LDA,*)
      REAL*8, INTENT(OUT) :: B(LDB,*)
      ! local variables
#ifndef _MKL_
      INTEGER :: NBLKSZ, LROWS, LCOLS, MAXROW, MAXCOL
      INTEGER :: I, J, IB, JB
#endif

      IF (NROWS.LE.0.OR.NCOLS.LE.0) THEN
        WRITE(5,'(1X,A)') 'DTRANS: Error: invalid dimension(s)'
        WRITE(5,'(1X,2(A,I9))') 'NROWS = ', NROWS, 'NCOLS = ', NCOLS
        CALL AbEnd
      ELSE IF(NROWS > LDA .OR. NCOLS > LDB) THEN
        WRITE(5,'(1X,A)') 'DTRANS: Error: dimension(s) out-of-bounds'
        WRITE(5,'(1X,2(A,I9))') 'NROWS = ', NROWS, 'NCOLS = ', NCOLS
        WRITE(5,'(1X,2(A,I9))') 'LDA   = ', LDA  , 'LDB   = ', LDB
        CALL AbEnd
      END IF

#ifdef _MKL_
      CALL mkl_domatcopy('C', 'T', NROWS, NCOLS, 1.0D0, A, LDA, B, LDB)
#else
      NBLKSZ=8
      LROWS=MOD(NROWS,NBLKSZ)
      LCOLS=MOD(NCOLS,NBLKSZ)
      MAXROW=NROWS-LROWS
      MAXCOL=NCOLS-LCOLS
      IF (MAXROW.GT.0.AND.MAXCOL.GT.0) THEN
        DO IB=1,MAXROW,NBLKSZ
          DO JB=1,MAXCOL,NBLKSZ
            DO I=IB,IB+NBLKSZ-1
              DO J=JB,JB+NBLKSZ-1
                B(J,I)=A(I,J)
              END DO
            END DO
          END DO
        END DO
      END IF
* remainder of the blocks
      IF (MAXROW.GT.0.AND.LCOLS.GT.0) THEN
        DO IB=1,MAXROW,NBLKSZ
          DO I=IB,IB+NBLKSZ-1
            DO J=MAXCOL+1,MAXCOL+LCOLS
              B(J,I)=A(I,J)
            END DO
          END DO
        END DO
      END IF
      IF (MAXCOL.GT.0.AND.LROWS.GT.0) THEN
        DO JB=1,MAXCOL,NBLKSZ
          DO I=MAXROW+1,MAXROW+LROWS
            DO J=JB,JB+NBLKSZ-1
              B(J,I)=A(I,J)
            END DO
          END DO
        END DO
      END IF
* remainder block
      IF (LROWS.GT.0.AND.LCOLS.GT.0) THEN
        DO I=MAXROW+1,MAXROW+LROWS
          DO J=MAXCOL+1,MAXCOL+LCOLS
            B(J,I)=A(I,J)
          END DO
        END DO
      END IF
#endif
      END SUBROUTINE
