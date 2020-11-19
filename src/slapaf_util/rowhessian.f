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
* Copyright (C) Giovanni Ghigo                                         *
************************************************************************
      Subroutine RowHessian(nIter,nInter,nRowH,mRowH,Delta,g)
      Implicit Real*8 (A-H,O-Z)
      Real*8 g(nInter,nIter)
      Integer mRowH(10)
#include "stdalloc.fh"
      Real*8 rDum(1)
      Real*8, Allocatable:: Hss_Q(:,:)
*
      Call mma_allocate(Hss_Q,nInter,nInter,Label='Hss_Q')
      Call Get_dArray('Hss_Q',Hss_Q,nInter**2)
      Call Put_dArray('Hss_upd',rDum,0)
*
      Call RowHessian_internal(nIter,nInter,nRowH,mRowH,Delta,Hss_Q,g)
      Write (6,*)
      Write (6,*) ' Partial numerical differentiation is finished!'
*
      Call Put_dArray('Hss_Q',Hss_Q,nInter**2)
      Call mma_deallocate(Hss_Q)
*
      Return
      End
      Subroutine RowHessian_Internal(nIter,nInter,nRowH,mRowH,Delta,H,g)
************************************************************************
*                                                                      *
* Object: Numerical estimation of single rows and columns of Hessian   *
* Called from: RlxCtl when lRowH=.True. & iter.EQ.NmIter               *
* Author: Giovanni Ghigo, University of Torino, Italy                  *
*                                                                      *
************************************************************************
      Implicit Real*8 (A-H,O-Z)
      Real*8 g(nInter,nIter),H(nInter,nInter)
      Integer mRowH(10)
#include "real.fh"
*
#ifdef _DEBUGPRINT_
      Write(6,*) 'RowHessian:'
      Call RecPrt('Initial Hessian',' ',H,nInter,nInter)
      Call RecPrt('Gradient      g:','(10F9.6)', g,nInter,nIter)
#endif
*
* --- Evaluate the Hessian
*
      Do iRowH = 1, nRowH
         iInter = mRowH(iRowH)
         Do jInter = 1, nInter
            H(iInter,jInter) = (g(jInter,1) - g(jInter,iRowH+1)) / Delta
            H(jInter,iInter) = H(iInter,jInter)
         End Do
      EndDo
*
* --- Symmetrize
*
      Do iInter = 1, nInter
         Do jInter = 1, nInter
            dElement = (H(iInter,jInter)+H(jInter,iInter))/Two
            H(iInter,jInter) = dElement
            H(jInter,iInter) = dElement
         EndDo
      EndDo
#ifdef _DEBUGPRINT_
      Call RecPrt('Final Hessian',' ',H,nInter,nInter)
#endif
      Return
      End
