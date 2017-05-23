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
      SUBROUTINE VDIV(A,LA,B,LB,C,LC,N)
      Implicit Real*8 (a-h,o-z)
      DIMENSION A(*),B(*),C(*)
      DO 10 I=0,N-1
         C(1+LC*I)=B(1+LB*I)/A(1+LA*I)
10    CONTINUE
      RETURN
      END
