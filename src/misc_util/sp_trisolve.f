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
*  Sp_TriSolve
*
*> @ingroup Sparse
*> @brief
*>   Solves a triangular linear system, with a sparse matrix
*> @author Ignacio Fdez. Galv&aacute;n
*>
*> @details
*> Solves the linear system \f$ A x = b \f$, where \p A is a sparse triangular matrix.
*> The side argument can be either ``'L'`` if \p A is lower triangular or ``'U'`` if
*> it is upper triangular.
*> On output the vector \p x contains the solution.
*>
*> @param[in]  n    Size of the system
*> @param[in]  side Type of system
*> @param[in]  A    Matrix in sparse format
*> @param[in]  ija  Index vector of matrix \p A
*> @param[in]  b    Vector of independent terms
*> @param[out] x    Solution vector
************************************************************************
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
