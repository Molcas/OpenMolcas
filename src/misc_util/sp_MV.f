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
*  Sp_MV
*
*> @ingroup Sparse
*> @brief
*>   Compute a matrix-vector product \f$ y \leftarrow \alpha A x + \beta y \f$, with a sparse matrix
*> @author Ignacio Fdez. Galv&aacute;n
*>
*> @details
*> Equivalent to ::DGeMV or ::DSyMV, with a sparse matrix \p A.
*>   \f[ y \leftarrow \alpha A x + \beta y \f]
*>
*> @param[in]     n     Size of the system
*> @param[in]     alpha Factor for the multiplication
*> @param[in]     A     Input matrix in sparse format
*> @param[in]     ija   Index vector of matrix \p A
*> @param[in]     x     Vector to multiply
*> @param[in]     beta  Factor for the initial vector
*> @param[in,out] y     Result vector
************************************************************************
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
