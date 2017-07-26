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
*  Allocate_Work
*
*> @brief
*>   Allocate real memory of length \p n and return the corresponding pointer \p ip
*> @author R. Lindh
*>
*> @details
*> Short parameter list wrapper to ::GETMEM for allocation of
*> memory of type ``REAL*8``. Label defaulted to '``AW``'.
*>
*> @param[out] ip pointer to \c Work
*> @param[in]  n  length
************************************************************************
      Subroutine Allocate_Work(ip,n)
      Implicit Real*8 (a-h,o-z)
#include "WrkSpc.fh"
*
      Call GetMem('AW','Allo','Real',ip,n)
*
      Return
      End
