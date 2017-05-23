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
      SUBROUTINE ORDER_AXIS(A,N)
      Implicit Real*8 (a-h,o-z)
      REAL*8 A(N)
      DO 10 I = 1, N-1
         DO 20 J = I+1, N
            IF (A(I).GT.A(J)) THEN
               SAVE=A(I)
               A(I)=A(J)
               A(J)=SAVE
            END IF
 20      CONTINUE
 10   CONTINUE
      RETURN
      END
