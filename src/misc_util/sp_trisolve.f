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
* Copyright (C) 2014, Ignacio Fdez. Galvan                             *
************************************************************************
*-----------------------------------------------------------------------
* <DOC>
*   <NAME>Sp\_TriSolve</NAME>
*   <SYNTAX>Call Sp\_TriSolve(n,side,A,ija,b,x)</Syntax>
*   <ARGUMENTS>
*     \Argument{n}{Size of the system}{Integer}{in}
*     \Argument{side}{Type of system}{Character}{in}
*     \Argument{A}{Matrix in sparse format}{Real*8 (*)}{in}
*     \Argument{ija}{Index vector of matrix A}{Integer (*)}{in}
*     \Argument{b}{Vector of independent terms}{Real*8 (n)}{in}
*     \Argument{x}{Solution vector}{Real*8 (n)}{out}
*   </ARGUMENTS>
*   <PURPOSE>Solves a triangular linear system, with a sparse matrix</PURPOSE>
*   <DEPENDENCIES></DEPENDENCIES>
*   <AUTHOR>I. Fdez. Galvan</AUTHOR>
*   <MODIFIED_BY></MODIFIED_BY>
*   <SIDE_EFFECTS></SIDE_EFFECTS>
*   <DESCRIPTION>
*     Solves the linear system $A x = b$, where A is a sparse triangular matrix.
*     The side argument can be either `L' if A is lower triangular or `U' if
*     it is upper triangular.
*     On output the vector x contains the solution.
*   </DESCRIPTION>
* </DOC>
*-----------------------------------------------------------------------
      SUBROUTINE Sp_TriSolve(n,side,A,ija,b,x)
      IMPLICIT NONE
      INTEGER n, ija(*), i, j, k
      REAL*8 A(*), b(n), x(n)
      CHARACTER side

      IF (side.EQ.'L') THEN
        DO i=1,n
          x(i)=b(i)
          DO k=ija(i),ija(i+1)-1
            j=ija(k)
            x(i)=x(i)-A(k)*x(j)
          END DO
          x(i)=x(i)/A(i)
        END DO
      ELSE IF (side.EQ.'U') THEN
        DO i=n,1,-1
          x(i)=b(i)
          DO k=ija(i),ija(i+1)-1
            j=ija(k)
            x(i)=x(i)-A(k)*x(j)
          END DO
          x(i)=x(i)/A(i)
        END DO
      END IF

      END
