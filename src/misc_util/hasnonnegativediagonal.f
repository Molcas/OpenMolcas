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
      Logical Function hasNonnegativeDiagonal(M,n)
************************************************************
*
*   <DOC>
*     <Name>hasNonnegativeDiagonal</Name>
*     <Syntax>hasNonnegativeDiagonal(M,n)</Syntax>
*     <Arguments>
*       \Argument{M}{n-by-n quadratic matrix}{Real*8 array}{in}
*       \Argument{n}{Dimension of M}{Integer}{in}
*     </Arguments>
*     <Purpose>Return .True. if all diagonal elements of M
*     are non-negative.
*     </Purpose>
*     <Dependencies></Dependencies>
*     <Author>Thomas Bondo Pedersen</Author>
*     <Modified_by></Modified_by>
*     <Side_Effects></Side_Effects>
*     <Description>Returns .True. if for all i, M(i,i) .ge. 0.0d0
*     </Description>
*    </DOC>
*
************************************************************
      Implicit None
      Integer n
      Real*8  M(n,n)

      Integer i

      hasNonnegativeDiagonal=.True.
      Do i=1,n
         If (M(i,i).lt.0.0d0) Then
            hasNonnegativeDiagonal=.False.
            Return
         End If
      End Do

      End
