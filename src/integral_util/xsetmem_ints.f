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
      subroutine xsetmem_ints(mem)
      implicit real*8 (a-h,o-z)
#include "status.fh"
c
      If (XMem_Status.eq.Active) Then
         Call WarningMessage(2,
     &               'External handling of scratch already active!')
         Call Abend()
      End If
      Call GetMem('SewMem','ALLO','REAL',ipt,mem)
      Call SetMem_Ints(ipt,Mem)
      XMem_Status=Active
*
      Return
      End
