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
      REAL*8 FUNCTION THETA(M,N)
      IMPLICIT REAL*8 (A-H,O-Z)
#include "crelop.fh"
C
C     INTEGRATION OVER THETA. INCLUDES A FACTOR SIN(TH)
C     FOR THE VOLUME ELEMENT
C
      IF (MOD(N,2).EQ.1) GOTO 10
      THETA=GA(M+2)*GA(N+1)/GA(M+N+3)
      RETURN
10    THETA=0.D0
      RETURN
      END
