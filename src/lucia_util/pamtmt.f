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
* Copyright (C) 1988, Jeppe Olsen                                      *
************************************************************************
      SUBROUTINE PAMTMT(X,T,SCR,NORB)
*
* GENERATE PER AKE'S T MATRIX FROM AN
* ORBITAL ROTATION MATRIX X
*
* T IS OBTAINED AS A STRICTLY LOWER TRIANGULAR
* MATRIX TL AND AN UPPER TRIANGULAR MATRIX TU
*
*         TL = 1 - L
*         TU = U ** -1
*
* WHERE L AND U ARISES FROM A LU DECOMPOSITION OF
* X :
*         X = L * U
* WITH L BEING A LOWER TRIANGULAR MATRIX WITH UNIT ON THE
* DIAGONAL AND U IS AN UPPER TRIANGULAR MATRIX
*
* JEPPE OLSEN OCTOBER 1988
*
      IMPLICIT REAL*8 (A-H,O-Z)
      Integer nOrb
      Real*8 X(NORB,NORB),T(NORB,NORB)
      Real*8 SCR(nOrb**2+nOrb*(nOrb+1)/2)
*
      NTEST = 00
      IF(NTEST.GE.2) THEN
        WRITE(6,*) ' Wellcome to PAMTMT '
        WRITE(6,*) ' =================='
        WRITE(6,*)
      END IF
*. Allocate local memory
      KLFREE = 1
C     KLL = KFLREE
      KLL = KLFREE
      KLFREE = KLL + NORB*(NORB+1)/2
      KLU = KLFREE
      KLFREE = KLU + NORB ** 2
*.LU factorize X
      CALL LULU(X,SCR(KLL),SCR(KLU),NORB)
*.Expand U to full matrix
      CALL SETVEC(T,0.0D0,NORB ** 2 )
      DO 10 I = 1,NORB
      DO 11 J = I,NORB
        T(I,J) = SCR(KLU-1+J*(J-1)/2+I)
   11 CONTINUE
   10 CONTINUE
      IF ( NTEST .GE. 100 ) THEN
        WRITE(6,*) ' MATRIX TO BE INVERTED '
        CALL WRTMAT(T,NORB,NORB,NORB,NORB)
      END IF
*.Invert U
      CALL INVMAT(T,SCR(KLU),NORB,NORB,ISING)
      IF ( NTEST .GE. 100 ) THEN
        WRITE(6,*) ' INVERTED MATRIX '
        CALL WRTMAT(T,NORB,NORB,NORB,NORB)
      END IF
*.Subtract L
      DO 20 I = 1, NORB
      DO 21 J = 1,I-1
       T(I,J)= - SCR(KLL-1+I*(I-1)/2+J)
   21 CONTINUE
   20 CONTINUE
*
      IF( NTEST .GE. 2 ) THEN
        WRITE(6,*) ' INPUT X MATRIX '
        CALL WRTMAT(X,NORB,NORB,NORB,NORB)
        WRITE(6,*) ' T MATRIX '
        CALL WRTMAT(T,NORB,NORB,NORB,NORB)
      END IF
*
      RETURN
      END
