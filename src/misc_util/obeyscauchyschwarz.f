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
      Logical Function obeysCauchySchwarz(M,n,Tol)
************************************************************
*
*   <DOC>
*     <Name>obeysCauchySchwarz</Name>
*     <Syntax>obeysCauchySchwarz(M,n,Tol)</Syntax>
*     <Arguments>
*       \Argument{M}{n-by-n quadratic matrix to test}{Real*8 array}{in}
*       \Argument{n}{Dimension of M}{Integer}{in}
*       \Argument{Tol}{Tolerance}{Real*8 scalar}{in}
*     </Arguments>
*     <Purpose>Return .True. if M obeys the Cauchy-Schwarz inequality
*              within tolerance Tol</Purpose>
*     <Dependencies></Dependencies>
*     <Author>Thomas Bondo Pedersen</Author>
*     <Modified_by></Modified_by>
*     <Side_Effects></Side_Effects>
*     <Description> Test if matrix M obeys the C-S inequality. If for
*     some i.gt.j,
*     abs(M(i,j)**2-M(i,i)*M(j,j)) .gt. Tol
*     return .False.
*     </Description>
*    </DOC>
*
************************************************************
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
