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
#include "stdalloc.fh"

      Integer i
      Integer, Allocatable:: LB(:)

      Call mma_allocate(LB,[0,nProcs-1],Label='LB')
      LB(:)=0

      LB(myRank) = nnBstRT(1)
      Call Cho_GAIGop(LB,nProcs,'+')
      Call Cho_Head('Cholesky vector dimension on each node','=',80,
     &              LuPri)
      Do i = 0,nProcs-1
         Write(LuPri,'(2X,A,I4,5X,A,I7)')
     &   'Node:',i,'Dimension:',LB(i)
      End Do

      Call mma_deallocate(LB)

      End
