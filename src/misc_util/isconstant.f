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
*  isConstant
*
*> @brief
*>   Return ``.True.`` if all elements of \p M are identical to \p Const within tolerance \p Tol
*> @author Thomas Bondo Pedersen
*>
*> @details
*> Test if array \p M is constant. If \f$ \forall i \; |M_i - \text{Const}| \le \text{Tol} \f$
*> return ``.True.``.
*>
*> @param[in] M     Array to test
*> @param[in] n     Dimension of \p M
*> @param[in] Const Constant
*> @param[in] Tol   Tolerance
*>
*> @return ``.True.`` if all elements of \p M are identical to \p Const within tolerance \p Tol
************************************************************************
      Logical Function isConstant(M,n,Const,Tol)
      Implicit None
      Integer n
      Real*8  M(n)
      Real*8  Const
      Real*8  Tol

      Integer i

      isConstant=.True.
      i=0
      Do While (i.lt.n .and. isConstant)
         i=i+1
         isConstant=abs(M(i)-Const).le.Tol
      End Do

      End
