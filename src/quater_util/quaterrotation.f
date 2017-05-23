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
      subroutine quaterRotation(Q,U,V)
************************************************************
*
*   <DOC>
*     <Name>quaterRotation</Name>
*     <Syntax>quaterRotation(Q,U,V)</Syntax>
*     <Arguments>
*       \Argument{Q}{Quaternion used for the rotation, Dimension(4)}{Real*8}{in}
*       \Argument{U}{Vector to be rotated, Dimension(3)}{Real*8}{in}
*       \Argument{V}{Vector rotated, Dimension(3)}{Real*8}{out}
*     </Arguments>
*     <Purpose>Performs the rotation of U in V via the quaternion Q</Purpose>
*     <Dependencies>quater util, blas and util</Dependencies>
*     <Author>Y. Carissan</Author>
*     <Modified_by></Modified_by>
*     <Side_Effects>none</Side_Effects>
*     <Description>
*        Performs the rotation of U in V via the quaternion Q
*     </Description>
*    </DOC>
*
************************************************************
      Implicit none
#include "WrkSpc.fh"
#include "debug.fh"
#include "real.fh"
      Real*8 U(3),V(3),Q(0:3)
      Real*8 T(3)
      Real*8 C1,C2,C3
      real*8 ddot_

      Call CheckQuater(Q)
      Call Cross(Q(1),U,T)    ! T=QxU
c
      C1=Two*Q(0)**2-One          ! C1 = 2 * Q(0)^2 - 1
      C2=2*Q(0)                   ! C2 = 2 * Q(0)
      C3=Two*ddot_(3,Q(1),1,U,1)   ! C3 = 2 * Q.U
      V(1) = C1*U(1) - C2*T(1) + C3*Q(1)
      V(2) = C1*U(2) - C2*T(2) + C3*Q(2)
      V(3) = C1*U(3) - C2*T(3) + C3*Q(3)
      Return
      End
