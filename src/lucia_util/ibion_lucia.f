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
      INTEGER FUNCTION IBION_LUCIA(M,N)
C PAM05:
C The following code does not always work.
C Replaced by call to my 'NOVERM' routine, renamed IBINOM.
C
C BIONOMIAL COEFFICIENT (M / N ) = IFAC(M)/(IFAC(M-N)*IFAC(N))
C
*      IB = 1
*      IF(M-N.GE.N) THEN
*         DO 100 K = (M-N+1), M
*           IB = IB * K
*  100    CONTINUE
*         IB = IB/IFAC(N)
*      ELSE
*         DO 200 K = N+1,M
*           IB = IB * K
*  200    CONTINUE
*         IB = IB/IFAC(M-N)
*      END IF
**
*      IBION_LUCIA = IB
*C
      IBION_LUCIA=IBINOM(M,N)
      RETURN
      END
