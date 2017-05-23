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
      SUBROUTINE AUXC(N,NBC,TC,ANSC)
      IMPLICIT REAL*8 (A-H,O-Z)
      X=1D0/(1D0+TC)
      Y=TC*X
      X=sqrt(X**(NBC+1))
      DUM=1D0
      IF (N.GT.1) THEN
        DO 10 I=N,2,-1
          DUM=DUM*Y*DBLE(NBC+2*I-3)/DBLE(2*I-2) + 1.D0
10      CONTINUE
      ENDIF
      ANSC=X*DUM
      RETURN
      END
