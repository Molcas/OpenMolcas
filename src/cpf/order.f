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
* Copyright (C) 1986, Per E. M. Siegbahn                               *
*               1986, Margareta R. A. Blomberg                         *
************************************************************************
      SUBROUTINE ORDER(C,D,N)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION C(*),D(*)
      IF(N.EQ.1)RETURN
      N1=N-1
      DO 10 I=1,N1
      I1=I+1
      DO 20 J=I1,N
      IF(D(I).LE.D(J))GO TO 20
      DT=D(I)
      D(I)=D(J)
      D(J)=DT
      IN=(I-1)*N
      IOUT=(J-1)*N
      DO 30 K=1,N
      CT=C(IN+K)
      C(IN+K)=C(IOUT+K)
      C(IOUT+K)=CT
30    CONTINUE
20    CONTINUE
10    CONTINUE
      RETURN
      END
