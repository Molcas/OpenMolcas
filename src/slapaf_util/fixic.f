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
      Subroutine FIXIC(nFix,SS,mInt,B,NDIM,F,Label,u)
      Implicit Real*8 (a-h,o-z)
#include "real.fh"
#include "WrkSpc.fh"
      Real*8 SS(mInt), B(nDim*mInt), F(nDim), u(nDim)
      Character*8 Label(mInt)
*
      Call qEnter('Fixic')
*
*     write out the internal coordinates which will be fixed
*
      WRITE (6,*)
      WRITE (6,*)
     &      ' Following internal coordinates are fixed'
      WRITE (6,*)
*
*     loop over all internal coordinates to be fixed
*
      Do 716 I = mInt-nFix+1,mInt
         WRITE (6,'(A,A,E10.3,A)')
     &             Label(i),' with a gradient of ',SS(I),
     &             ' is frozen and the gradient is annihilated'
         SS(i) = Zero
716   Continue
*
*     now transform remaining internal coordinates back to cartesian ba
*                          -1 +
*                    fx = u  B  fq
*
      Call GetMem('uInv','Allo','Real',ipuInv,nDim**2)
      call dcopy_(nDim**2, [Zero],0,Work(ipuInv),1)
      Do i = 1, nDim
         ii=(i-1)*nDim + i
          Work(ipuInv+ii-1)=One/u(i)
      End Do
*     Call RecPrt('uInv',' ',Work(ipuInv),nDim,nDim)
      Call GetMem('uB','Allo','Real',ipuB,mInt*nDim)
      Call DGEMM_('N','N',
     &            nDim,mInt,nDim,
     &            1.0d0,Work(ipuInv),nDim,
     &            B,nDim,
     &            0.0d0,Work(ipuB),nDim)
*     Call RecPrt('uInvB',' ',Work(ipuB),nDim,mInt)
      Call DGEMM_('N','N',
     &            nDim,1,mInt,
     &            1.0d0,Work(ipuB),nDim,
     &            SS,mInt,
     &            0.0d0,F,nDim)
*     Call RecPrt('F',' ',F,mInt,1)
      Call GetMem('uB','Free','Real',ipuB,mInt*nDim)
      Call GetMem('uInv','Free','Real',ipuInv,nDim**2)
*
      Call qExit('Fixic')
      Return
      End
