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
      REAL*8 FUNCTION SIMPLM(NMP,F,X)
C
C...  Calculates  Intg/xmin,xmax/(F(x) dx)
C     if the F function has been evaluated in a logarithmic
C     mesh.
C
C     Note that the integral equals Ingt/xmin,xmax/(F.x d(ln x)),
C     which is what is actually calculated.
C
C
C     NMP:    number of mesh points
C     F:      array with the function values
C     X:      array with the variable values of the log. mesh
C
C
      IMPLICIT REAL*8 (A-H,O-Z)
      Real*8    F(NMP),X(NMP)
C
C...  checks the precision of the incr. of the log.mesh
      DLM= LOG(X(2))-LOG(X(1))
      DO 11 J=2,5
        DELTA= LOG(X(J+1))-LOG(X(J))
        IF (abs(DELTA-DLM).LT.1D-8) GO TO 11
          WRITE (6,602)
          Call Abend
11    CONTINUE
C
      IF (MOD(NMP,2).NE.0) THEN
        N=NMP
      ELSE
        N=NMP-1
      ENDIF
C...  Odd number of points
      SUM=0D0
      NM2=N-2
      DO 10 I=1,NM2,2
        IP1=I+1
        IP2=I+2
        SUM=SUM + F(I)*X(I)+4D0*F(IP1)*X(IP1)+F(IP2)*X(IP2)
10    CONTINUE
      SIMPLM=SUM*DLM/3D0
      IF (N.EQ.NMP) RETURN
C
C...  even number of points: add up the remaining
      SUM=2.5D0*F(NMP)*X(NMP)+4D0*F(N)*X(N)-0.5D0*F(N-1)*X(N-2)
      SIMPLM=SIMPLM + SUM*DLM/6D0
      RETURN
602   FORMAT (' SIMPLM: Increment of the log mesh not constant')
      END
