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
*  getMatrix
*
*> @brief
*>   Utility routine for quaternion resolution
*> @author Y. Carissan
*>
*> @details
*> Fills the output matrix with the rotation matrix.
*>
*> @param[out] M Rotation matrix at output
************************************************************************
      subroutine getMatrix(M)
      Implicit none
#include "WrkSpc.fh"
#include "rotation.fh"
      Real*8 M(3,3)

      if (matrixSet) then
        call dcopy_(6,RotMatrix,1,M,1)
      else
        Call SysAbendMsg("RdInput","Rotation Matrix was not set",
     &          "")
      end if

      Return
      End
