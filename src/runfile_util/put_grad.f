************************************************************************
* This file is part of OpenMolcas.                                     *
*                                                                      *
* OpenMolcas is free software; you can redistribute it and/or modify   *
* it under the terms of the GNU Lesser General Public License, v. 2.1. *
* OpenMolcas is distributed in the hope that it will be useful, but it *
* is provided "as is" and without any express or implied warranties.   *
* For more details see the full text of the license in the file        *
* LICENSE or in <http://www.gnu.org/licenses/>.                        *
************************************************************************
      Subroutine Put_Grad(Grad,nGrad)
      Implicit Real*8 (a-h,o-z)
#include "SysDef.fh"

      Real*8       Grad(nGrad)
      Character*24 Label

      Label='GRAD'
      Call Put_dArray(Label,Grad,nGrad)
*
      Call Get_iScalar('Grad ready',iGO)
      iGO = iOr(iGO,2**0)
      Call Put_iScalar('Grad ready',iGO)

      Return
      End
