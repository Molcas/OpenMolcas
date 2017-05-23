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
      SubRoutine Cho_P_iSySh(ip_iSySh)
C
C     Purpose: allocate integer array iSySh(nnShl_G)
C
C     (placed in separate routine as global data is accessed)
C
      Implicit None
      Integer ip_iSySh
#include "choglob.fh"

      Integer l_iSySh

      l_iSySh = nnShl_G
      Call GetMem('DiaSh','Allo','Real',ip_iSySh,l_iSySh)

      End
