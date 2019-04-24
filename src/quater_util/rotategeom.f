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
*  rotategeom
*
*> @brief
*>   Performs the rotation of \p Geom1 with the \p Q quaternion and stores the result in \p Geom2
*> @author Y. Carissan
*>
*> @details
*> Performs the rotation of \p Geom1 with the \p Q quaternion
*> and stores the result in \p Geom2.
*>
*> @param[in]  Q     Input Quaternion
*> @param[in]  nat   Number of atoms
*> @param[in]  Geom1 Geometry to be rotated, xyz coordinates
*> @param[out] Geom2 Output geometry, xyz coordinates
************************************************************************
      subroutine rotategeom(Q, nat, Geom1, Geom2)
      implicit none
      Integer iat,nat
      Real*8 Q(0:3)
      Real*8 V(3)
      Real*8 Geom1(3,nat),Geom2(3,nat)

      Do iat=1,nat
        call dcopy_(3,Geom1(1,iat),1,V,1)
        Call QuaterRotation(Q,V,Geom2(1,iat))
      End do
      Return
      End
