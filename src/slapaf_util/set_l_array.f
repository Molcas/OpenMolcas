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
      Subroutine set_l_Array(Array_l,nInter,BaseLine,Hessian)
      Implicit Real*8 (a-h,o-z)
      Real*8 Array_l(nInter), Hessian(nInter,nInter), l_max
*
      Call RecPrt('set_l_Array: Hessian',' ',Hessian,nInter,nInter)
      Write (6,*) 'BaseLine=',BaseLine
*
*     Gives a Kriging Hessian for a single point of Kriging with
*     a diagonal which is identical to the diagonal values of
*     the HMF ad hoc Hessian.
*
*     l_max=80.0D0
*     l_max=100.0D0
      Do i = 1, nInter
*
*        Make sure that the characteristic length is not too long.
*
         Hss=Max(Abs(Hessian(i,i)),0.050D0)
         Hss=Max(Abs(Hessian(i,i)),0.010D0)
         Hss=Max(Abs(Hessian(i,i)),0.020D0)
         Hss=Max(Abs(Hessian(i,i)),0.025D0)
         Array_l(i)=Sqrt((5.0D0*BaseLine)/(3.0D0*Hss))
*        Hss=Abs(Hessian(i,i))
*        Array_l(i)=Min(l_max,Sqrt((5.0D0*BaseLine)/(3.0D0*Hss)))
*
      End Do
      Call RecPrt('Array_l',' ',Array_l,1,nInter)
*
      Return
      End
