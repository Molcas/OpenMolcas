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
      SubRoutine Cho_P_DiaSh(ip_DiaSh)
C
C     Purpose: allocate double precision array DiaSh(nnShl_G)
C
C     (placed in separate routine as global data is accessed)
C
      Implicit None
      Integer ip_DiaSh
#include "choglob.fh"

      Integer l_DiaSh

      l_DiaSh = nnShl_G
      Call GetMem('DiaSh','Allo','Real',ip_DiaSh,l_DiaSh)

      End
