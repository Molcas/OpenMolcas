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
* Copyright (C) 1994, Roland Lindh                                     *
************************************************************************
      Subroutine QNR(nInter,nIter,dq,H,g)
************************************************************************
*                                                                      *
*     Author: Roland Lindh, Dept. of Theoretical Chemistry,            *
*             University of Lund, SWEDEN                               *
*             December '94                                             *
************************************************************************
      Implicit Real*8 (A-H,O-Z)
#include "print.fh"
#include "real.fh"
#include "WrkSpc.fh"
      Real*8 dq(nInter,nIter), H(nInter,nInter), g(nInter,nIter+1)
*
*     Call QEnter('QNR')
      iRout = 115
      iPrint = nPrint(iRout)
*
*-----Compute a new independent geometry by relaxation of
*     the gradient vector.
*
!     Call Allocate_iWork(ipTmp,1)
      Call Allocate_iWork(ipTmp,nInter*nInter)
!     iWork(ipTmp)=0
      Call IZero(iWork(ipTmp),nInter*nInter)
!     Call DZero(nInter,dq(1,nIter))
      Call DZero(dq(1,nIter),nInter)
      Call CG_Solver(nInter,nInter**2,H,iWork(ipTmp),
     &               g(1,nIter),dq(1,nIter),Info)
      Call Free_iWork(ipTmp)
      If (Info.lt.0) Then
        Call SysAbendMsg('QNR','Conjugate gradients not converged',' ')
      End If
*
      If (iPrint.ge.99) Then
        Write(6,*) 'CG converged in ',Info,' iterations'
        Call RecPrt(' H ',' ',H,nInter,nInter)
        Call RecPrt(' dq',' ',dq(1,nIter),1,nInter)
      End If
*
*     Call QExit('QNR')
      Return
      End
