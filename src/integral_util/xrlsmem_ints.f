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
      Subroutine xRlsMem_Ints
      implicit real*8 (a-h,o-z)
#include "status.fh"
c
      If (XMem_Status.eq.InActive) Then
C        Write (6,*)
C        Write (6,*) 'No external scratch handling to deactivate!'
C        Write (6,*)
         Return
      End If
      Call RlsMem_Ints
      XMem_Status=InActive
*
      Return
      End
