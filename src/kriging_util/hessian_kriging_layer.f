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
      Subroutine Hessian_Kriging_Layer(qInt,Hessian,nInter)
      Use Limbo
      Implicit None
#include "stdalloc.fh"
      Integer nInter
      Real*8 qInt(nInter), Hessian(nInter,nInter)
      Real*8, Allocatable:: qInt_s(:), Hessian_s(:,:)
*#define _DEBUGPRINT_
#ifdef _DEBUGPRINT_
      Call RecPrt('KHL: qInt',' ',qInt,1,nInter)
#endif
*
      Call mma_allocate(qInt_s,nInter,Label='qInt_s')
      Call mma_allocate(Hessian_s,nInter,nInter,Label='Hessian_s')
*
      Call Trans_K(U,qInt,qInt_s,nInter,1)
      Call Hessian_kriging(qInt_s,Hessian,nInter)
      Call BackTrans_K (U,Hessian,Hessian_s,nInter,nInter)
      Call BackTrans_K2(U,Hessian_s,Hessian,nInter,nInter)
*
      Call mma_deallocate(Hessian_s)
      Call mma_deallocate(qInt_s)
#ifdef _DEBUGPRINT_
      Call RecPrt('KHL: Hessian',' ',Hessian,nInter,nInter)
#endif
*
      End Subroutine Hessian_Kriging_Layer
