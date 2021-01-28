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
      SubRoutine Cho_SubScr_Final()
C
C     Purpose: finalize (de-allocate memory) screening in vector
C              subtraction.
C
      use ChoSubScr, only: DSubScr, ip_DSPNm, l_DSPNm
      Implicit None
#include "stdalloc.fh"

      If (l_DSPNm .gt. 0) Then
         Call Cho_Mem('DSPNm','Free','Real',ip_DSPNm,l_DSPNm)
         l_DSPNm = 0
      End If

      If (Allocated(DSubScr)) Call mma_deallocate(DSubScr)

      End
