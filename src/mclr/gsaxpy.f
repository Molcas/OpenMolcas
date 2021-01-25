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
* Copyright (C) 1994, Jeppe Olsen                                      *
************************************************************************
      SUBROUTINE GSAXPY(AB,A,B,NABCOL,NACOL,NROW,IABCOL,IACOL)
*
* AB(I,IABCOL(J)) = AB(I,IABCOL(J)) + A(I,IACOL(K))*B(K,J)
*
*
* Jeppe Olsen, Spring of 94 Daughter of MSAXPY*
*

      IMPLICIT REAL*8 (A-H,O-Z)
*. Input
      DIMENSION A(NROW,*),B(NACOL,NABCOL)
      DIMENSION IACOL(NACOL),IABCOL(NABCOL)
*. Output
      DIMENSION AB(NROW,*)
*
      IWAY = 2
C      ICRAY = 0
      IF(IWAY.EQ.1) THEN
*. Straightforward sequence of SAXPY's
         DO 1100 J = 1, NABCOL
           DO 1200 K = 1, NACOL
             JACT = IABCOL(J)
             KACT = IACOL(K)
             FACTOR = B(K,J)
C             IF(ICRAY.EQ.1) THEN
C               CALL SAXPY(NROW,FACTOR,A(1,KACT),1,AB(1,JACT),1)
C             ELSE
               DO 1000 I = 1, NROW
                 AB(I,JACT) = AB(I,JACT) + FACTOR*A(I,KACT)
 1000          CONTINUE
C             END IF
 1200      CONTINUE
 1100    CONTINUE
*
      ELSE IF (IWAY .EQ. 2 ) THEN
*. Unrolling over columns of A
       NROL = 5
       NRES = MOD(NACOL,NROL)
*. overhead
       IF(NRES.EQ.1) THEN
         DO 2201 J = 1, NABCOL
          JACT = IABCOL(J)
          K1ACT = IACOL(1)
          B1J = B(1,J)
C          IF(ICRAY.EQ.1) THEN
C            CALL SAXPY(NROW,B1J,A(1,K1ACT),1,AB(1,JACT),1)
C          ELSE
            DO 2101 I = 1, NROW
              AB(I,JACT) = AB(I,JACT) + A(I,K1ACT)*B1J
 2101       CONTINUE
C          END IF
 2201    CONTINUE
       ELSE IF (NRES .EQ. 2 ) THEN
         DO 2202 J = 1, NABCOL
           K1ACT = IACOL(1)
           K2ACT = IACOL(2)
           B1J   = B(1,J)
           B2J   = B(2,J)
           JACT =  IABCOL(J)
           DO 2102 I = 1, NROW
            AB(I,JACT) = AB(I,JACT)
     &    + A(I,K1ACT)*B1J + A(I,K2ACT)*B2J
 2102      CONTINUE
 2202    CONTINUE
       ELSE IF (NRES .EQ. 3 ) THEN
         DO 2203 J = 1, NABCOL
           K1ACT = IACOL(1)
           K2ACT = IACOL(2)
           K3ACT = IACOL(3)
           JACT =  IABCOL(J)
           B1J   = B(1,J)
           B2J   = B(2,J)
           B3J   = B(3,J)
           DO 2103 I = 1, NROW
            AB(I,JACT) = AB(I,JACT)
     &    + A(I,K1ACT)*B1J + A(I,K2ACT)*B2J
     &    + A(I,K3ACT)*B3J
 2103      CONTINUE
 2203    CONTINUE
       ELSE IF (NRES .EQ. 4 ) THEN
         DO 2204 J = 1, NABCOL
           K1ACT = IACOL(1)
           K2ACT = IACOL(2)
           K3ACT = IACOL(3)
           K4ACT = IACOL(4)
           JACT =  IABCOL(J)
           B1J   = B(1,J)
           B2J   = B(2,J)
           B3J   = B(3,J)
           B4J   = B(4,J)
           DO 2104 I = 1, NROW
            AB(I,JACT) = AB(I,JACT)
     &    + A(I,K1ACT)*B1J + A(I,K2ACT)*B2J
     &    + A(I,K3ACT)*B3J + A(I,K4ACT)*B4J
 2104      CONTINUE
 2204    CONTINUE
        END IF
*. ( End of Overhead )
        DO 2305 K = NRES+1,NACOL,NROL
          DO 2205 J = 1, NABCOL
            K1ACT = IACOL(K)
            K2ACT = IACOL(K+1)
            K3ACT = IACOL(K+2)
            K4ACT = IACOL(K+3)
            K5ACT = IACOL(K+4)
            JACT =  IABCOL(J)
            B1J   = B(K,J)
            B2J   = B(K+1,J)
            B3J   = B(K+2,J)
            B4J   = B(K+3,J)
            B5J   = B(K+4,J)
            DO 2105 I = 1, NROW
             AB(I,JACT) = AB(I,JACT)
     &     + A(I,K1ACT)*B1J + A(I,K2ACT)*B2J
     &     + A(I,K3ACT)*B3J + A(I,K4ACT)*B4J
     &     + A(I,K5ACT)*B5J
 2105       CONTINUE
 2205     CONTINUE
 2305   CONTINUE
      END IF
*( End of IWAY branching )
*
      RETURN
      END
