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
* Copyright (C) 2019, Roland Lindh                                     *
************************************************************************
      Real*8 Function Restriction_Dispersion(q,dq,nInter)
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
      Implicit Real*8 (a-h,o-z)
#include "stdalloc.fh"
      Integer nInter
      Real*8 q(nInter), dq(nInter)
      Real*8, Allocatable:: qNext(:)
*
*define _DEBUG_
      Call mma_Allocate(qNext,nInter,Label='qNext')
      Do i = 1, nInter
         qNext(i)=q(i)+dq(i)
      End Do
#ifdef _DEBUG_
      Call RecPrt('q',' ',q,nInter,1)
      Call RecPrt('dq',' ',dq,nInter,1)
      Call RecPrt('qNext',' ',qNext,nInter,1)
#endif
      Call Dispersion_Kriging(qNext,y,nInter)
      Restriction_Dispersion=y
      Call mma_Deallocate(qNext)
*
      Return
      End
