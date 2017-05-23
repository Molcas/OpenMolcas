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
      SUBROUTINE REOR_MAT(Y,X,I1,I2,J1,J2)
      implicit none
      integer I1,I2,J1,J2, I,J,IJ,L
      REAL*8 Y,X,ONE,SCALAR
      PARAMETER (ONE=1.D0)
      DIMENSION Y(I1,I2,J1,J2),X(*)
C
CSUBROUTINE EX24 EXCHANGE OF THE 2nd and 4th index in a matrix
C
      ENTRY EX24(Y,X,I1,I2,J1,J2)
      IJ=1
      DO J=1,I2
         DO L=1,J1
            DO I=1,J2
               CALL DCOPY_(I1,Y(1,J,L,I),1,X(IJ),1)
               IJ=IJ+I1
            ENDDO
         ENDDO
      ENDDO
      RETURN
C
C
      ENTRY EX24_A(Y,X,I1,I2,J1,J2,SCALAR)
C
CS EX24_A EXCHANGE OF THE 2nd and 4th index in a matrix+ADD
C
      IJ=1
      DO J=1,I2
         DO L=1,J1
            DO I=1,J2
               CALL DAXPY_(I1,SCALAR,Y(1,J,L,I),1,X(IJ),1)
!              CALL DCOPY_(I1,Y(1,J,L,I),1,X(IJ),1)
               IJ=IJ+I1
            ENDDO
         ENDDO
      ENDDO
      RETURN
C
C
      ENTRY EX34_A(Y,X,I1,I2,J1,J2,SCALAR)
C
CS EX34_A EXCHANGE OF THE 3nd and 4th index in a matrix+ADD
C
      IJ=1
!     DO J=1,I2
      DO L=1,J1
         DO I=1,J2
            CALL DAXPY_(I1*I2,SCALAR,Y(1,1,L,I),1,X(IJ),1)
!           CALL DCOPY_(I1,Y(1,J,L,I),1,X(IJ),1)
            IJ=IJ+I1*I2
         ENDDO
      ENDDO
!     ENDDO
      RETURN
C
C
      ENTRY EX23(Y,X,I1,I2,J1,J2)
C
      IJ=1
      DO J=1,J2
         DO L=1,I2
            DO I=1,J1
               CALL DCOPY_(I1,Y(1,L,I,J),1,X(IJ),1)
               IJ=IJ+I1
            ENDDO
         ENDDO
      ENDDO
      RETURN
C
C
      ENTRY EX23_A(Y,X,I1,I2,J1,J2)
C
      IJ=1
      DO J=1,J2
         DO L=1,I2
            DO I=1,J1
               CALL DAXPY_(I1,ONE,Y(1,L,I,J),1,X(IJ),1)
               IJ=IJ+I1
            ENDDO
         ENDDO
      ENDDO
      RETURN
C
C
      ENTRY EX312(Y,X,I1,I2,J1,J2)
C
      IJ=1
      DO J=1,J2
         DO L=1,I2
            DO I=1,I1
               CALL DCOPY_(J1,Y(I,L,1,J),I1*I2,X(IJ),1)
               IJ=IJ+J1
            ENDDO
         ENDDO
      ENDDO
      RETURN
C
C
      ENTRY EX2413(Y,X,I1,I2,J1,J2)
C
      IJ=1
      DO J=1,J1
         DO L=1,I1
            DO I=1,J2
               CALL DCOPY_(I2,Y(L,1,J,I),I1,X(IJ),1)
               IJ=IJ+I2
            ENDDO
         ENDDO
      ENDDO

      RETURN
C
C
      ENTRY EX2413_A(Y,X,I1,I2,J1,J2)
C
      IJ=1
      DO J=1,J1
         DO L=1,I1
            DO I=1,J2
               CALL DAXPY_(I2,ONE,Y(L,1,J,I),I1,X(IJ),1)
               IJ=IJ+I2
            ENDDO
         ENDDO
      ENDDO
      RETURN
C
C
      ENTRY EX423(Y,X,I1,I2,J1,J2)
C
      IJ=1
      DO J=1,J1
         DO L=1,I2
            DO I=1,J2
               CALL DCOPY_(I1,Y(1,L,J,I),1,X(IJ),1)
               IJ=IJ+I1
            ENDDO
         ENDDO
      ENDDO
      RETURN
C
C
      ENTRY EX423_A(Y,X,I1,I2,J1,J2,SCALAR)
C
      IJ=1
      DO J=1,J1
         DO L=1,I2
            DO I=1,J2
               CALL DAXPY_(I1,SCALAR,Y(1,L,J,I),1,X(IJ),1)
C              CALL DCOPY_(I1,Y(1,L,J,I),1,X(IJ),1)
               IJ=IJ+I1
            ENDDO
         ENDDO
      ENDDO
      RETURN
      END
