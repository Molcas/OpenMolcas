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
      Subroutine QuaterSetup(U1,U2,V1,V2)
************************************************************
*
*   <DOC>
*     <Name>QuaterSetup</Name>
*     <Syntax>Call QuaterSetup(U1,U2,V1,V2)</Syntax>
*     <Arguments>
*       \Argument{U1}{Input Vector, Dimension(3)}{Real*8}{in}
*       \Argument{V1}{Input Vector, Dimension(3)}{Real*8}{in}
*       \Argument{U2}{Input Vector, Dimension(3)}{Real*8}{in}
*       \Argument{V2}{Input and output Vector, Dimension(3)}{Real*8}{in-out}
*     </Arguments>
*     <Purpose>
*        Transforms V2 such that U1.V1 = U2.V2
*       </Purpose>
*     <Dependencies></Dependencies>
*     <Author>Y. Carissan</Author>
*     <Modified_by></Modified_by>
*     <Side_Effects>none</Side_Effects>
*     <Description>
*              Transforms V2 such that U1.V1 = U2.V2
*     </Description>
*    </DOC>
*
************************************************************
      Implicit none
#include "debug.fh"
#include "real.fh"
      Real*8 U1(3),V1(3)
      Real*8 U2(3),V2(3)
      Real*8 coeff,U1dU2,V1dV2
      real*8 ddot_
      external ddot_
*
      call normalizeVec(U1)
      call normalizeVec(V1)
      call normalizeVec(U2)
      call normalizeVec(V2)
*
      if (debug) then
        call RecPrt("IN QUATERSETUP normalized U1","",U1,3,1)
        call RecPrt("IN QUATERSETUP normalized V1","",V1,3,1)
        call RecPrt("IN QUATERSETUP normalized U2","",U2,3,1)
        call RecPrt("IN QUATERSETUP normalized V2","",V2,3,1)
      end if
*
      U1dU2=ddot_(3,U1,1,U2,1)
      V1dV2=ddot_(3,V1,1,V2,1)
      coeff = ( One-U1dU2**2 )/( One-V1dV2**2 )
      coeff = sqrt(coeff)
      V2(1) = ( U1dU2 - V1dV2 * coeff ) * V1(1) + coeff * V2(1)
      V2(2) = ( U1dU2 - V1dV2 * coeff ) * V1(2) + coeff * V2(2)
      V2(3) = ( U1dU2 - V1dV2 * coeff ) * V1(3) + coeff * V2(3)
      if (debug) call RecPrt("IN QUATERSETUP modified V2","",V2,3,1)
      Return
      End
