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
      SubRoutine Sq2Tri(Sq,Tri,n)
************************************************************
*
*   <DOC>
*     <Name>Sq2Tri</Name>
*     <Syntax>Call Sq2Tri(Sq,Tri,n)</Syntax>
*     <Arguments>
*       \Argument{Sq}{Square storage array, dim. Sq(n,n)}{Real*8}{in}
*       \Argument{Tri}{Lower triangular storage array, dim.
*                      Tri(n*(n+1)/2}{Real*8}{out}
*       \Argument{n}{Dimension}{Integer}{in}
*     </Arguments>
*     <Purpose>Convert from square to lower triangular storage</Purpose>
*     <Dependencies></Dependencies>
*     <Author>Thomas Bondo Pedersen</Author>
*     <Modified_by></Modified_by>
*     <Side_Effects></Side_Effects>
*     <Description>
*        Perform the extraction
*           Tri(i*(i-1)/2+j) = Sq(i,j)
*        where i >= j.
*     </Description>
*    </DOC>
*
************************************************************
      Implicit None
      Integer n
      Real*8  Sq(n,n), Tri(n*(n+1)/2)

      Integer i,j, iTri

      iTri(i,j)=i*(i-1)/2 + j

      Do j = 1,n
         Do i = j,n
            Tri(iTri(i,j)) = Sq(i,j)
         End Do
      End Do

      End
