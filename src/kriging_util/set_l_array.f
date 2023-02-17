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
      Subroutine set_l_Array(Array_l,nInter,BaseLine,Hessian,HDiag)
      Implicit Real*8 (a-h,o-z)
      Integer nInter
      Real*8 Array_l(nInter)
      Real*8, Optional:: Hessian(nInter,nInter)
      Real*8, Optional:: HDiag(nInter)
      Real*8 H_min
*
!#define _DEBUGPRINT_
#ifdef _DEBUGPRINT_
      Call RecPrt('set_l_Array: Hessian',' ',Hessian,nInter,nInter)
      Write (6,*) 'BaseLine=',BaseLine
#endif
*
*     Gives a Kriging Hessian for a single point of Kriging with
*     a diagonal which is identical to the diagonal values of
*     the HMF ad hoc Hessian.
*
      H_min=0.025D0
      If (Present(Hessian)) Then
         Do i = 1, nInter
*
*           Make sure that the characteristic length is not too long.
*
            Hss=Max(Abs(Hessian(i,i)),H_min)
            Hessian(i,i)=Hss
            Array_l(i)=Sqrt((5.0D0*BaseLine)/(3.0D0*Hss))
*
         End Do
      Else
         Do i = 1, nInter
*
*           Make sure that the characteristic length is not too long.
*
            Hss=Max(Abs(HDiag(i)),H_min)
            HDiag(i)=Hss
            Array_l(i)=Sqrt((5.0D0*BaseLine)/(3.0D0*Hss))
*
         End Do
      End If
*
#ifdef _DEBUGPRINT_
      Write (6,*) 'H_min=',H_Min
      Call RecPrt('Array_l',' ',Array_l,1,nInter)
#endif
*
      Return
      End
