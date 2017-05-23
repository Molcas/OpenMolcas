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
      SubRoutine Cho_PTS_Stat()
      Implicit None
#include "cholesky.fh"
#include "choptr.fh"
#include "choglob.fh"
#include "cho_para_info.fh"
#include "WrkSpc.fh"

      Integer iTmp

      If (l_IntMap.lt.1) Then
         l_IntMap=nnShl
         Call GetMem('IntMap','Allo','Inte',ip_IntMap,l_IntMap)
         iTmp=0
         Call IDAFile(LuMap,2,iWork(ip_IntMap),l_IntMap,iTmp)
      End If

      If (Cho_Real_Par) Then
         Call iSwap(nSym,NumCho,1,NumCho_G,1)
         iTmp = NumChT
         NumChT = NumChT_G
         Call Cho_Stat()
         NumChT = iTmp
         Call iSwap(nSym,NumCho,1,NumCho_G,1)
      Else
         Call Cho_Stat()
      End If

      If (l_IntMap.gt.0) Then
         Call GetMem('IntMap','Free','Inte',ip_IntMap,l_IntMap)
      End If

      End
