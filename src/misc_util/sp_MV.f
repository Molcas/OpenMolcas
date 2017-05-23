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
*   <NAME>Sp\_MV</NAME>
*   <SYNTAX>Call Sp\_MV(n,alpha,A,ija,x,beta,y)</Syntax>
*   <ARGUMENTS>
*     \Argument{n}{Size of the system}{Integer}{in}
*     \Argument{alpha}{Factor for the multiplication}{Real*8}{in}
*     \Argument{A}{Input matrix in sparse format}{Real*8 (*)}{in}
*     \Argument{ija}{Index vector of matrix A}{Integer (*)}{in}
*     \Argument{x}{Vector to multiply}{Real*8 (n)}{in}
*     \Argument{beta}{Factor for the initial vector}{Real*8}{in}
*     \Argument{y}{Result vector}{Real*8 (n)}{inout}
*   </ARGUMENTS>
*   <PURPOSE>Compute a matrix-vector product y = alpha*A*x + beta*y, with a sparse matrix</PURPOSE>
*   <DEPENDENCIES></DEPENDENCIES>
*   <AUTHOR>I. Fdez. Galvan</AUTHOR>
*   <MODIFIED_BY></MODIFIED_BY>
*   <SIDE_EFFECTS></SIDE_EFFECTS>
*   <DESCRIPTION>
*     Equivalent to dgemv or dsymv, with a sparse matrix A.
*       $y = \alpha A x + \beta y$
*   </DESCRIPTION>
* </DOC>
*-----------------------------------------------------------------------
      SUBROUTINE Sp_MV(n,alpha,A,ija,x,beta,y)
      IMPLICIT NONE
      INTEGER n, ija(*), i, j, k
      REAL*8 alpha, beta, A(*), x(n), y(n)
#include "real.fh"

c
c     Very simple routine, but split in different cases
c     to gain efficiency
      IF (A(n+1).GT.Zero) THEN
        IF ((beta.EQ.Zero).AND.(alpha.EQ.One)) THEN
          DO i=1,n
            y(i)=A(i)*x(i)
            DO k=ija(i),ija(i+1)-1
              j=ija(k)
              y(i)=y(i)+A(k)*x(j)
              y(j)=y(j)+A(k)*x(i)
            END DO
          END DO
        ELSE
          DO i=1,n
            y(i)=beta*y(i)+alpha*A(i)*x(i)
            DO k=ija(i),ija(i+1)-1
              j=ija(k)
              y(i)=y(i)+alpha*A(k)*x(j)
              y(j)=y(j)+alpha*A(k)*x(i)
            END DO
          END DO
        END IF
      ELSE
        IF ((beta.EQ.Zero).AND.(alpha.EQ.One)) THEN
          DO i=1,n
            y(i)=A(i)*x(i)
            DO k=ija(i),ija(i+1)-1
              y(i)=y(i)+A(k)*x(ija(k))
            END DO
          END DO
        ELSE
          DO i=1,n
            y(i)=beta*y(i)+alpha*A(i)*x(i)
            DO k=ija(i),ija(i+1)-1
              y(i)=y(i)+alpha*A(k)*x(ija(k))
            END DO
          END DO
        END IF
      END IF

      END
