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
      subroutine rotategeom(Q, nat, Geom1, Geom2)
************************************************************
*
*   <DOC>
*     <Name>rotategeom</Name>
*     <Syntax>Call rotategeom(Q, nat, Geom1, Geom2)</Syntax>
*     <Arguments>
*       \Argument{Q}{Input Quaternion, Dimension(4)}{Real*8}{in}
*       \Argument{nat}{Number of atoms}{Integer}{in}
*       \Argument{Geom1}{Geometry to be rotated, xyz coordinates, Dimension(nat,3)}{Real*8}{in}
*       \Argument{Geom2}{Output geometry, xyz coordinates, Dimension(nat,3)}{Real*8}{out}
*     </Arguments>
*     <Purpose>Performs the rotation of Geom1 with the Q quaternion
*       and stores the result in Geom2</Purpose>
*     <Dependencies>quater util and blas</Dependencies>
*     <Author>Y. Carissan</Author>
*     <Modified_by></Modified_by>
*     <Side_Effects>none</Side_Effects>
*     <Description>
*     Performs the rotation of Geom1 with the Q quaternion
*       and stores the result in Geom2
*     </Description>
*    </DOC>
*
************************************************************

      implicit none
      Integer iat,nat
      Real*8 Q(0:4)
      Real*8 V(3)
      Real*8 Geom1(3,nat),Geom2(3,nat)

      Do iat=1,nat
        call dcopy_(3,Geom1(1,iat),1,V,1)
        Call QuaterRotation(Q,V,Geom2(1,iat))
      End do
      Return
      End
