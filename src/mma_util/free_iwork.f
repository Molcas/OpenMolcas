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
* Copyright (C) Roland Lindh                                           *
************************************************************************
*  Free_iWork
*
*> @brief
*>   Deallocate memory in \c iWork associated with pointer \p ip
*> @author Roland Lindh
*>
*> @details
*> The array associated to pointer/identifier \p ip in \c iWork is deallocated.
*>
*> @param[in,out] ip pointer to memory in \c iWork
************************************************************************
      Subroutine Free_iWork(ip)
      Implicit Real*8 (a-h,o-z)
#include "WrkSpc.fh"
*
      Call GetMem('AiW','Free','Inte',ip,nDum)
*
      Return
      End
