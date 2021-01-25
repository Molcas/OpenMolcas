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
      SubRoutine Cho_PrintLB()
      Use Para_Info, Only: MyRank, nProcs
      Implicit None
#include "cholesky.fh"
#include "WrkSpc.fh"

      Integer ip_LB, l_LB, i

      l_LB = nProcs
      Call GetMem('LoadB','Allo','Inte',ip_LB,l_LB)

      Call iZero(iWork(ip_LB),l_LB)
      iWork(ip_LB+myRank) = nnBstRT(1)
      Call Cho_GAIGop(iWork(ip_LB),l_LB,'+')
      Call Cho_Head('Cholesky vector dimension on each node','=',80,
     &              LuPri)
      Do i = 0,nProcs-1
         Write(LuPri,'(2X,A,I4,5X,A,I7)')
     &   'Node:',i,'Dimension:',iWork(ip_LB+i)
      End Do

      Call GetMem('LoadB','Free','Inte',ip_LB,l_LB)

      End
