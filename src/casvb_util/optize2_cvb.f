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
      subroutine optize2_cvb(fx,nparm,ioptc,
     >  dx,grad,iter_is_1,
     >  opta,optb)
      implicit real*8 (a-h,o-z)
      logical opth,skipupd,first_time
      logical iter_is_1,close2conv_begin
      logical close2conv,converged,wrongstat,scalesmall1
#include "opt_cvb.fh"
#include "locopt1_cvb.fh"
#include "locopt2_cvb.fh"
#include "trst_cvb.fh"
#include "tune_cvb.fh"

#include "formats_cvb.fh"
      external opta,optb
      dimension dx(nparm),grad(nparm)
      save zero,close2conv,converged
      data zero/0.d0/

      converged=.false.
      if(iter_is_1)close2conv=.false.

c  << Initialization >>
      call grad_cvb(grad)
      call ddproj_cvb(grad,nparm)
      grdnrm=dnrm2_(nparm,grad,1)
      call opta(nparm)

      if(ip.ge.2)write(6,'(/a)')' *****   2. order optimizer   *****'

c  << Now trust region control >>
      exp_tc=exp
      first_time=.true.
      opth=.false.
      iopth=0
100   call trust_cvb(iopth,opth,maxize,fx,fxbest,exp,
     >  hh,dxnrm,ioptc,scalesmall1,close2conv,converged,skipupd)
      if(ioptc.eq.-2)return

c    << Make update >>
      if(.not.(skipupd.or.hh.eq.zero))then

200     close2conv_begin=close2conv

        call optb(nparm,dxnrm,grdnrm,close2conv)

        if(first_time)then
          first_time=.false.

          call testconv_cvb(fx,nparm,
     >      dx,grad,exp_tc,
     >      close2conv,converged,wrongstat)

          if(close2conv.and..not.close2conv_begin)goto 200
        endif
        if((ip.eq.2.and..not.opth).or.ip.ge.3)then
          s11=ddot_(nparm,dx,1,dx,1)
          s22=ddot_(nparm,grad,1,grad,1)
          s12=ddot_(nparm,dx,1,grad,1)
          write(6,formAD)
     >      ' Overlap between normalized vectors <DX|GRAD> :',
     >      s12/sqrt(s11*s22)
        endif
        call fxdx_cvb(fx,.false.,dx)
      endif
      if(opth)goto 100

      if(ioptc.ge.-1.and.hh.ne.zero)then
        if(ip.ge.2)then
          write(6,'(a)')' '
          write(6,formAF)' HH & norm of update :',hh,dxnrm
        endif
        call update_cvb(dx)
      endif
      if(converged)then
        ioptc=0
      elseif(close2conv.and.endifclose)then
        ioptc=-3
      else
        ioptc=1
      endif
      return
      end
