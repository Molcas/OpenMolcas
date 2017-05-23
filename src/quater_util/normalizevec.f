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
* Copyright (C) Yannick Carissan                                       *
************************************************************************
      subroutine normalizeVec(V)
************************************************************
*
*   <DOC>
*     <Name>NormalizeVec</Name>
*     <Syntax>Call NormalizeVec(V)</Syntax>
*     <Arguments>
*       \Argument{V}{Vector to be normalized}{Real*8}{in-out}
*     </Arguments>
*     <Purpose>Normalize the vector V</Purpose>
*     <Dependencies>blas</Dependencies>
*     <Author>Y. Carissan</Author>
*     <Modified_by></Modified_by>
*     <Side_Effects>none</Side_Effects>
*     <Description>
*        Normalize the vector V
*     </Description>
*    </DOC>
*
************************************************************

      Implicit none
#include "WrkSpc.fh"
#include "debug.fh"

      Real*8 V(3)
      Real*8 norm,dnrm2_
      External dnrm2_
      Integer i

      norm=dnrm2_(3,V,1)
      Do i=1,3
        V(i)=V(i)/norm
      End Do

      Return
      End
