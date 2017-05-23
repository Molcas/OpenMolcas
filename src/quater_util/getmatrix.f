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
      subroutine getMatrix(M)
************************************************************
*
*   <DOC>
*     <Name>getMatrix(M)</Name>
*     <Syntax>getMatrix(M)</Syntax>
*     <Arguments>
*       \Argument{M}{Rotation matrix at output, Dimension(3,3)}{Real*8}{out}
*     </Arguments>
*     <Purpose>Utility routine for quaternion resolution</Purpose>
*     <Dependencies>blas setmatrix must have been called</Dependencies>
*     <Author>Y. Carissan</Author>
*     <Modified_by></Modified_by>
*     <Side_Effects>none</Side_Effects>
*     <Description>
*        Fills the output matrix with the rotation matrix.
*     </Description>
*    </DOC>
*
************************************************************
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
