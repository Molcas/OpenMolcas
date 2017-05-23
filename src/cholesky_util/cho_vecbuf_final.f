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
      SubRoutine Cho_VecBuf_Final()
C
C     Purpose: deallocate and finalize vector buffer.
C
      Implicit None
#include "cholesky.fh"
#include "chovecbuf.fh"

      If (l_ChVBuf .gt. 0) Then
         Call Cho_Mem('CHVBUF','Free','Real',ip_ChVBuf,l_ChVBuf)
      End If
      l_ChVBuf=0
      ip_ChVBuf=0
      If (l_ChVBfI .gt. 0) Then
         Call Cho_Mem('ChVBfI','Free','Real',ip_ChVBfI,l_ChVBfI)
      End If
      l_ChVBfI=0
      ip_ChVBfI=0
      Call Cho_iZero(ip_ChVBuf_Sym,nSym)
      Call Cho_iZero(l_ChVBuf_Sym,nSym)
      Call Cho_iZero(ip_ChVBFI_Sym,nSym)
      Call Cho_iZero(l_ChVBFI_Sym,nSym)
      Call Cho_iZero(nVec_in_Buf,nSym)

      End
