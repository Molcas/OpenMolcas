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
      SUBROUTINE SQUAR2(A,N)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION A(N,N)
      DO 10 I=1,N
      NI=N-I+1
      CALL DCOPY_(NI,A(I,I),1,A(I,I),N)
10    CONTINUE
      RETURN
      END
