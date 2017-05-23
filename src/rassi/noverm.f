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
* Copyright (C) 1998, Per Ake Malmqvist                                *
************************************************************************
      INTEGER FUNCTION NOVERM(N,M)
      IMPLICIT REAL*8 (A-H,O-Z)
      SAVE INIT,NOMTAB
      DIMENSION NOMTAB(225)
      DATA INIT / 0 /
      NOVERM=0
      IF(N.LT.0) RETURN
      MM=M
      IF(2*MM.GT.N) MM=N-M
      IF(MM.LT.0) RETURN
      NOVERM=1
      IF(MM.EQ.0) RETURN
      NOVERM=N
      IF(MM.EQ.1) RETURN
      IF(INIT.EQ.0) THEN
        IPOS=0
        DO I=4,32
         X=DBLE(I)
         DO J=2,I/2
          IPOS=IPOS+1
          X=(X*DBLE(I+1-J))/DBLE(J)
          NOMTAB(IPOS)=NINT(X)
         END DO
        END DO
CPAM      write(*,'(1x,5I9)') NOMTAB
        INIT=1
      END IF
      IF(N.LE.32) THEN
        NOVERM=NOMTAB(((N-3)**2)/4+MM-1)
      ELSE
        X=DBLE(NOVERM)
        DO K=2,MM
         X=(X*DBLE(N+1-K))/DBLE(K)
        END DO
        NOVERM=NINT(X)
        IF(X.NE.DBLE(NOVERM)) THEN
          WRITE(6,*)' NOVERM: Unable to compute N over M'
          WRITE(6,*)' N=',N
          WRITE(6,*)' M=',M
          CALL ABEND()
        END IF
      END IF
      RETURN
      END
