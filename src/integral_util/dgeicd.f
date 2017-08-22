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
      SUBROUTINE DGEICD(A,LDA,N,IOPT,RCOND,DET,AUX,NAUX)
C
C     COMPUTE THE INVERSE, THE RECIPROCAL CONDITION NUMBER AND THE
C     DETERMINANT OF MATRIX A
C
      IMPLICIT REAL*8 (A-H,O-Z)
#include "WrkSpc.fh"
      DIMENSION A(LDA,N),AUX(NAUX),DET(2)
      DIMENSION IPVT(16),Z(16),TEMP2(2)
C
      IF ( N.LT.0 .OR. N.GT.LDA ) THEN
         WRITE(6,*)
         WRITE(6,'(6X,A)') '*** ERROR IN SUBROUTINE DGEICD ***'
         WRITE(6,'(6X,A)') 'ORDER OF MATRIX A IS OUT OF BOUNDS'
         WRITE(6,*)
         Call Abend()
      END IF
      IF ( IOPT.LT.0 .OR. IOPT.GT.3 ) THEN
         WRITE(6,*)
         WRITE(6,'(6X,A)') '*** ERROR IN SUBROUTINE DGEICD ***'
         WRITE(6,'(6X,A)') '    OPTION KEY IS OUT OF BOUNDS'
         WRITE(6,*)
         Call Abend()
      END IF
      IF ( N.GE.16 .AND. NAUX.LT.2*N ) THEN
         WRITE(6,*)
         WRITE(6,'(6X,A)') '*** ERROR IN SUBROUTINE DGEICD ***'
         WRITE(6,'(6X,A)') '      WORK AREA IS TO SMALL'
         WRITE(6,*)
         Call Abend()
      END IF
C
      IF ( N.LT.16 ) THEN
         CALL DGECO(A,LDA,N,IPVT,TEMP1,Z)
      ELSE
C        CALL DGECO(A,LDA,N,AUX(1),TEMP1,AUX(N+1))
*
*------- This trick to have the correct data type for the 4th argument.
*
         ipAux=ip_of_iWork(Aux(1))
         CALL DGECO(A,LDA,N,iWork(ipAux),TEMP1,AUX(N+1))
      END IF
      IF ( (1.0D0+TEMP1).EQ.1.0D0 ) THEN
         WRITE(6,*)
         WRITE(6,'(6X,A)') '*** ERROR IN SUBROUTINE DGEICD ***'
         WRITE(6,'(6X,A)') '      THIS A SINGULAR MATRIX'
         WRITE(6,*)
         Call Abend()
      END IF
      IF ( IOPT.EQ.1 .OR. IOPT.EQ.3 ) RCOND=TEMP1
      JOB=11
      IF ( N.LT.16 ) THEN
         CALL DGEDI(A,LDA,N,IPVT,TEMP2,Z,JOB)
      ELSE
C        CALL DGEDI(A,LDA,N,AUX(1),TEMP2,AUX(N+1),JOB)
*
*------- This trick to have the correct data type for the 4th argument.
*
         ipAux=ip_of_iWork(Aux(1))
         CALL DGEDI(A,LDA,N,iWork(ipAux),TEMP2,AUX(N+1),JOB)
      END IF
      IF ( IOPT.EQ.2 .OR. IOPT.EQ.3 ) THEN
         DET(1)=TEMP2(1)
         DET(2)=TEMP2(2)
      END IF
C
      RETURN
      END
