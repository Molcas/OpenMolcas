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
* Copyright (C) Thomas Bondo Pedersen                                  *
************************************************************************
*  obeysCauchySchwarz
*
*> @brief
*>   Return ``.True.`` if \p M obeys the Cauchy--Schwarz inequality within tolerance \p Tol
*> @author Thomas Bondo Pedersen
*>
*> @details
*> Test if matrix \p M obeys the Cauchy--Schwarz inequality.
*> If \f$ \exists i \gt j \; |M_{ij}^2 - M_{ii}M_{jj}| \gt \text{Tol} \f$
*> return ``.False.``.
*>
*> @param[in] M   \p n &times; \p n square matrix to test
*> @param[in] n   Dimension of \p M
*> @param[in] Tol Tolerance
*>
*> @return ``.True.`` if \p M obeys the Cauchy--Schwarz inequality within tolerance \p Tol
************************************************************************
      Logical Function obeysCauchySchwarz(M,n,Tol)
      Implicit None
      Integer n
      Real*8  M(n,n)
      Real*8  Tol

      Integer i, j

      obeysCauchySchwarz=.True.
      Do j=1,n
         Do i=j+1,n
            If (M(i,j)**2 .gt. M(i,i)*M(j,j) .and.
     &          abs(M(i,j)**2-M(i,i)*M(j,j)).gt.Tol) Then
               obeysCauchySchwarz=.False.
               Return
            End If
         End Do
      End Do

      End
