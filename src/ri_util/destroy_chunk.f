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
      Subroutine Destroy_Chunk(ip_Chunk,ip_iMap)
      Integer ip_Chunk
#include "para_info.fh"
*
#ifdef _MOLCAS_MPP_
      If (Is_Real_Par()) Then
         Call GA_Destroy(ip_Chunk)
         Call Free_iWork(ip_iMap)
      Else
         Call Free_Work(ip_Chunk)
      End If
#else
      Call Free_Work(ip_Chunk)
c Avoid unused argument warnings
      If (.False.) Call Unused_integer(ip_iMap)
#endif
      Return
      End
