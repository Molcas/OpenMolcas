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
*  Free_Work
*
*> @brief
*>   Deallocate memory in \c Work associated with pointer \p ip
*> @author Roland Lindh
*>
*> @details
*> The array associated to pointer/identifier \p ip in \c Work is deallocated.
*>
*> @param[in,out] ip pointer to memory in \c Work
************************************************************************
      Subroutine Free_Work(ip)
      Implicit Real*8 (a-h,o-z)
#include "WrkSpc.fh"
*
      Call GetMem('AW','Free','Real',ip,0)
*
      Return
      End
