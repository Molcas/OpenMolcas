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
      Subroutine set_l_Array(Array_l,nInter)
      Implicit Real*8 (a-h,o-z)
#include "real.fh"
#include "stdalloc.fh"
      Real*8 Array_l(nInter)
      Real*8, Allocatable:: Hessian(:,:)
*
      Call mma_Allocate(Hessian,nInter,nInter,Label='Hessian')
      Call Mk_Hss_Q()
      Call Get_dArray('Hss_Q',Hessian,nInter**2)
*     Call RecPrt('Hessian',' ',Hessian,nInter,nInter)
*                                                                      *
      rmax=Zero
      Do i = 1, nInter
         iCase=1
         If (iCase.eq.1) Then
            Array_l(i)=One/Sqrt(Abs(Hessian(i,i)))
         Else If (iCase.eq.2) Then
            Array_l(i)=One/Abs(Hessian(i,i))
         Else If (iCase.eq.3) Then
            Array_l(i)=Abs(Hessian(i,i))
         Else If (iCase.eq.4) Then
            Array_l(i)=Sqrt(Abs(Hessian(i,i)))
         Else
            Write (6,*) 'set_l_array: illegal iCase value'
            Call Abend()
         End If
         If (Array_l(i).gt.rmax) rmax=Array_l(i)
      End Do
*     Call RecPrt('Raw',' ',Array_l,1,nInter)
      Call DScal_(nInter,One/rmax,Array_l,1)
*     Call RecPrt('Scaled',' ',Array_l,1,nInter)
*
      Call mma_Deallocate(Hessian)
*
      Return
      End
