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
      Subroutine set_l_Array(Array_l,nInter,BaseLine)
      Implicit Real*8 (a-h,o-z)
#include "real.fh"
#include "stdalloc.fh"
      Real*8 Array_l(nInter)
      Real*8, Allocatable:: Hessian(:,:)
*
      Call mma_Allocate(Hessian,nInter,nInter,Label='Hessian')
      Call Mk_Hss_Q()
      Call Get_dArray('Hss_Q',Hessian,nInter**2)
*     Call RecPrt('set_l_Array: Hessian',' ',Hessian,nInter,nInter)
*
*     Gives a Kriging Hessian for a single point of Kriging with
*     a diagonal which is identical to the diagonal values of
*     the HMF ad hoc Hessian.
*
      Do i = 1, nInter
*
         Array_l(i)=Sqrt((5.0D0*BaseLine)/(3.0D0*Abs(Hessian(i,i))))
*
      End Do
*     Call RecPrt('Array_l',' ',Array_l,1,nInter)
*
      Call mma_Deallocate(Hessian)
*
      Return
      End
