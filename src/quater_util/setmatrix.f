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
*  setMatrix
*
*> @brief
*>   Computes the rotation matrix of the rotation stored in the input quaternion
*> @author Y. Carissan
*>
*> @details
*> Computes the rotation matrix of the rotation stored in the input quaternion.
*>
*> @param[in] Q Input Quaternion
************************************************************************
      subroutine setMatrix(Q)
      implicit real*8 (a-h,o-z)
#include "WrkSpc.fh"
#include "rotation.fh"
#include "debug.fh"
#include "real.fh"
      Real*8 Q(0:3)

      RotMatrix(1,1)=Two*(Q(0)*Q(0) + Q(1)*Q(1)) - One
      RotMatrix(1,2)=Two*(Q(1)*Q(2) + Q(0)*Q(3))
      RotMatrix(1,3)=Two*(Q(1)*Q(3) - Q(0)*Q(2))

      RotMatrix(2,1)=Two*(Q(1)*Q(2) - Q(0)*Q(3))
      RotMatrix(2,2)=Two*(Q(0)*Q(0) + Q(2)*Q(2)) - One
      RotMatrix(2,3)=Two*(Q(2)*Q(3) + Q(0)*Q(1))

      RotMatrix(3,1)=Two*(Q(1)*Q(3) + Q(0)*Q(2))
      RotMatrix(3,2)=Two*(Q(2)*Q(3) - Q(0)*Q(1))
      RotMatrix(3,3)=Two*(Q(0)*Q(0) + Q(3)*Q(3)) - One

      matrixSet=.true.
      if(debug) Call RecPrt("Rotation Matrix",' ',RotMatrix,3,3)

      Return
      End
