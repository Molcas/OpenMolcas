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
* Copyright (C) 2020, Roland Lindh                                     *
************************************************************************
      Real*8 Function Restriction_Disp_Con(x,dx,mInter)
************************************************************************
*                                                                      *
*     Object: External routine to evaluate a general constraint,       *
*             to be used in a constrained optimization. In this case   *
*             the constraint is a step size constraint.                *
*                                                                      *
*     Author: Roland Lindh, Dept. of Chemistry - BMC                   *
*             University of Uppsala, Sweden                            *
*                                                                      *
************************************************************************
      use RDC
      Implicit Real*8 (a-h,o-z)
#include "stdalloc.fh"
      Integer mInter
      Real*8 x(mInter), dx(mInter)
      Real*8, Allocatable:: qNext(:), du(:), dq(:)
*
      nInter=nInter_
      nLambda=nLambda_
      If (nInter-nLambda.ne.mInter) Then
         Write (6,*) 'RDC: Sanity check failed.'
         Call Abend()
      End If
      Call mma_Allocate(qNext,nInter,Label='qNext')
      qNext(:)=0.0D0
      Call mma_Allocate(du,nInter,Label='du')
      du(:)=0.0D0
      Call mma_Allocate(dq,nInter,Label='dq')
      dq(:)=0.0D0
*
*     Call RecPrt('dx',' ',dx,1,mInter)
      du(nLambda+1:mInter)=dx(1:mInter)
      du(1:nLambda)=dy_(1:nLambda)
*     Call RecPrt('du',' ',du,1,nInter)
*     Call RecPrt('T_',' ',T_,nInter,nInter)
      Call DGEMM_('N','N',
     &            nInter,1,nInter,
     &            1.0D0,T_,nInter,
     &            du,nInter,
     &            0.0D0,dq,nInter)
*     Call RecPrt('dq',' ',dq,1,nInter)
*     Call RecPrt('q_',' ',q_,1,nInter)
*
      Do i = 1, nInter
         qNext(i)=q_(i)+dq(i)
      End Do
*     Call RecPrt('qNext',' ',qNext,1,nInter)
*
      Call Dispersion_Kriging_Layer(qNext,y,nInter)
      Restriction_Disp_Con=y
*
      Call mma_Deallocate(dq)
      Call mma_Deallocate(du)
      Call mma_Deallocate(qNext)
*
      If (.False.) Then
         Call Unused_real_array(x)
      End If
      Return
      End
