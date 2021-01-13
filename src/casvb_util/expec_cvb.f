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
* Copyright (C) 1996-2006, Thorstein Thorsteinsson                     *
*               1996-2006, David L. Cooper                             *
************************************************************************
      subroutine expec_cvb(dxp,gradp,heigval,
     >  nnegeig,npr,exp,exp1,exp2)
      implicit real*8 (a-h,o-z)
      dimension dxp(npr),gradp(npr),heigval(npr)
      save zero,half
      data zero/0.d0/,half/0.5d0/

      exp1l=zero
      do 100 i=1,nnegeig
      exp1l=exp1l+dxp(i)*(gradp(i)+half*dxp(i)*heigval(i))
100   continue
      exp2l=zero
      do 200 i=nnegeig+1,npr
      exp2l=exp2l+dxp(i)*(gradp(i)+half*dxp(i)*heigval(i))
200   continue
      exp=exp1l+exp2l
      exp1=exp1l
      exp2=exp2l
      return
      end
      subroutine getdxp_cvb(dxp,gradp,heigval,nnegeig,npr,alfa)
      implicit real*8 (a-h,o-z)
      dimension dxp(npr),gradp(npr),heigval(npr)
      do 300 i=1,nnegeig
      dxp(i)=-gradp(i)/(heigval(i)-alfa)
300   continue
      do 400 i=nnegeig+1,npr
      dxp(i)=-gradp(i)/(heigval(i)+alfa)
400   continue
      return
      end
