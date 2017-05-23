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
      Logical Function isConstant(M,n,Const,Tol)
************************************************************
*
*   <DOC>
*     <Name>isConstant</Name>
*     <Syntax>isConstant(M,n,Const,Tol)</Syntax>
*     <Arguments>
*       \Argument{M}{Array to test}{Real*8 array}{in}
*       \Argument{n}{Dimension of M}{Integer}{in}
*       \Argument{Const}{Constant}{Real*8 scalar}{in}
*       \Argument{Tol}{Tolerance}{Real*8 scalar}{in}
*     </Arguments>
*     <Purpose>Return .True. if all elements of M are identical to Const
*              within tolerance Tol.
*     </Purpose>
*     <Dependencies></Dependencies>
*     <Author>Thomas Bondo Pedersen</Author>
*     <Modified_by></Modified_by>
*     <Side_Effects></Side_Effects>
*     <Description> Test if array M is constant. If for all i
*     abs(M(i)-Const) .le. Tol
*     return .True.
*     </Description>
*    </DOC>
*
************************************************************
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
