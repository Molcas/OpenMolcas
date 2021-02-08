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
* Copyright (C) 1994, Per Ake Malmqvist                                *
*               2012, Thomas Bondo Pedersen                            *
************************************************************************
      SUBROUTINE TRNSPS(N,M,A,B)
      IMPLICIT REAL*8 (A-H,O-Z)
C Return B = Transpose of A.
      DIMENSION A(N,M),B(M,N)
      DO 10 I=1,N
        DO 11 J=1,M
          B(J,I)=A(I,J)
  11    CONTINUE
  10  CONTINUE
      RETURN
      END
      Subroutine iTrnsps(n,m,a,b)
C Thomas Bondo Pedersen, August 2012: integer version of trnsps
      Implicit None
      Integer n
      Integer m
      Integer a(n,m)
      Integer b(m,n)
      Integer i, j
      Do i=1,n
         Do j=1,m
            b(j,i)=a(i,j)
         End Do
      End Do
      Return
      End
