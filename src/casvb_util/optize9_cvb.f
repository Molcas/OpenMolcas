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
* Copyright (C) 1996-2006, T. Thorsteinsson and D. L. Cooper           *
************************************************************************
      subroutine optize9_cvb(fx1,nparm,ioptc,
     >  hessdx,grad,dx)
      implicit real*8 (a-h,o-z)
#include "formats_cvb.fh"
      dimension hessdx(nparm),grad(nparm),dx(nparm)
      save zero,tenth,half,one
      data zero/0d0/,tenth/1d-1/,half/.5d0/,one/1d0/

      call grad_cvb(grad)

      dum=rand_cvb(.777d0)
      do 100 iparm=1,nparm
100   dx(iparm)=rand_cvb(zero)-half
      call nize_cvb(dx,1,dum,nparm,0,0)
      call fmove(dx,hessdx,nparm)
      call hess_cvb(hessdx)

      write(6,'(2a)')' Simple check of gradient and Hessian using ',
     >  'a random update vector :'
      e1=ddot_(nparm,dx,1,grad,1)
      e2=ddot_(nparm,dx,1,hessdx,1)
      write(6,'(a)')' '
      write(6,formChk1)' First-order change  :',e1
      write(6,formChk1)' Second-order change :',e2
      write(6,'(a)')' '

      write(6,formChk2)'Norm     ','DFX(act) ','DFX(pred)','Ratio    ',
     >  'F2(act)'
      cn=one
      do 200 it=1,10
      call fxdx_cvb(fx,.false.,dx)
      write(6,formChk3)cn,fx-fx1,cn*e1+cn*cn*half*e2,
     >  (fx-fx1)/(cn*e1+cn*cn*half*e2),(fx-fx1-cn*e1)/(cn*cn*half)
      call dscal_(nparm,tenth,dx,1)
200   cn=tenth*cn

      ioptc=0
      return
      end
