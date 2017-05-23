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
      Subroutine Get_lArray(Label,Log,nLog)
      Character*(*) Label
#include "WrkSpc.fh"
      Logical Log(nLog)
*
      Call Allocate_iWork(ip_Tmp,nLog)
      Call Get_iArray(Label,iWork(ip_Tmp),nLog)
      Call Int2Log(iWork(ip_Tmp),Log,nLog)
      Call Free_iWork(ip_Tmp)
*
      Return
      End
