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
      subroutine testconv2_cvb(close2conv,converged,wrongstat,
     >  act,zz,step,grad,npr,eigmn,eigmx,eigmna,nposeig,nnegeig)
      implicit real*8 (a-h,o-z)
#include "formats_cvb.fh"
      logical close2conv,converged,wrongstat,close_old
      logical dfx_is_small,step_is_small,grad_is_small
      logical correct_index,zz_ok
#include "tols_cvb.fh"
#include "print_cvb.fh"
      dimension step(npr),grad(npr)

      close_old=close2conv

      step_nrm=dnrm2_(npr,step,1)
      step_rms=step_nrm/sqrt(DBLE(npr))
      call findamx_cvb(step,npr,step_amx,idum)

      grad_nrm=dnrm2_(npr,grad,1)
      grad_rms=grad_nrm/sqrt(DBLE(npr))
      call findamx_cvb(grad,npr,grad_amx,idum)

      if(ip(3).ge.2)then
        if(nnegeig.gt.0.and.nposeig.gt.0)then
          write(6,*)' Maximum eigenvalue : ',eigmx,' of first ',
     >      nnegeig,' values.'
          write(6,*)' Minimum eigenvalue : ',eigmn,' of last  ',
     >      nposeig,' values.'
        elseif(nnegeig.gt.0)then
          write(6,*)' Maximum eigenvalue : ',eigmx
        else
          write(6,*)' Minimum eigenvalue : ',eigmn
        endif
      endif

c  Local region / final convergence reached if :
c  a) Change in F(X) small
c  b) Step is small
c  c) Gradient is small
c  d) Hessian has correct index
c  e) ACT/EXP ratio is within acceptable interval

      if(eigmna.gt.singul(1))then
        mm=1
      else
        mm=2
      endif
      dfx_is_small=(act.lt.dfx(mm))
      step_is_small=(step_nrm.lt.dx(2,mm).and.
     >  step_rms.lt.dx(3,mm).and.step_amx.lt.dx(1,mm))
      grad_is_small=(grad_nrm.lt.grd(2,mm).and.
     >  grad_rms.lt.grd(3,mm).and.grad_amx.lt.grd(1,mm))
      correct_index=eigmx.lt.sign(mm).and.eigmn.gt.-sign(mm)
      zz_ok=(zz.gt.zzmin(mm).and.zz.lt.zzmax(mm))

      close2conv=(dfx_is_small.and.step_is_small.and.grad_is_small
     >  .and.correct_index.and.zz_ok)

      if(eigmna.gt.singul(2))then
        mm=3
      else
        mm=4
      endif
      dfx_is_small=(act.lt.dfx(mm))
      step_is_small=(step_nrm.lt.dx(2,mm).and.
     >  step_rms.lt.dx(3,mm).and.step_amx.lt.dx(1,mm))
      grad_is_small=(grad_nrm.lt.grd(2,mm).and.
     >  grad_rms.lt.grd(3,mm).and.grad_amx.lt.grd(1,mm))
      correct_index=eigmx.lt.sign(mm).and.eigmn.gt.-sign(mm)
      zz_ok=(zz.gt.zzmin(mm).and.zz.lt.zzmax(mm))

      if(ip(3).ge.2)then
        write(6,'(/,a)')' Test of convergence :'
        write(6,'(a)')' ---------------------'
        call cvprt_cvb(' 1) Change in F(x) :',dfx_is_small)
        call cvprt_cvb(' 2) Step length    :',step_is_small)
        call cvprt_cvb(' 3) Grad norm      :',grad_is_small)
        call cvprt_cvb(' 4) Hessian index  :',correct_index)
        call cvprt_cvb(' 5) Act/Exp ratio  :',zz_ok)
        write(6,*)' '
        call cvprt2_cvb(' F(x) change   :',act,dfx(mm),1)
        call cvprt2_cvb(' Norm of step  :',step_nrm,dx(2,mm),1)
        call cvprt2_cvb(' RMS of step   :',step_rms,dx(3,mm),1)
        call cvprt2_cvb(' AMAX of step  :',step_amx,dx(1,mm),1)
        call cvprt2_cvb(' Norm of grad  :',grad_nrm,grd(2,mm),1)
        call cvprt2_cvb(' RMS of grad   :',grad_rms,grd(3,mm),1)
        call cvprt2_cvb(' AMAX of grad  :',grad_amx,grd(1,mm),1)
        call cvprt2_cvb(' Max. eigval   :',eigmx,sign(mm),1)
        call cvprt2_cvb(' Min. eigval   :',eigmn,-sign(mm),2)
        call cvprt2_cvb(' Act/Exp ratio :',zz,zzmin(mm),2)
        call cvprt2_cvb(' Act/Exp ratio :',zz,zzmax(mm),1)
        write(6,*)' '
      endif

      converged=(dfx_is_small.and.step_is_small.and.grad_is_small
     >  .and.correct_index.and.zz_ok)
      converged=(converged.and.close2conv)

      if(eigmna.gt.singul(3))then
        mm=5
      else
        mm=6
      endif
      dfx_is_small=(act.lt.dfx(mm))
      step_is_small=(step_nrm.lt.dx(2,mm).and.
     >  step_rms.lt.dx(3,mm).and.step_amx.lt.dx(1,mm))
      grad_is_small=(grad_nrm.lt.grd(2,mm).and.
     >  grad_rms.lt.grd(3,mm).and.grad_amx.lt.grd(1,mm))
      correct_index=eigmx.lt.sign(mm).and.eigmn.gt.-sign(mm)
      zz_ok=(zz.gt.zzmin(mm).and.zz.lt.zzmax(mm))

c  Wrong stationary point if otherwise converged, but Hessian
c  hasn't got correct index:
      wrongstat=(dfx_is_small.and.step_is_small.and.grad_is_small
     >  .and.(.not.correct_index).and.zz_ok)

      if(ip(3).ge.1.and.close2conv.and..not.(converged.or.close_old))
     >  write(6,'(a)')' Optimization entering local region.'
      if(converged.and.ip(3).ge.0)then
        write(6,formAD)
     >    ' Converged ... maximum update to coefficient:',step_amx
        if(eigmna.le.singul(2))then
          write(6,'(a)')' Warning - singular hessian!'
          write(6,formAD)' Smallest Hessian eigenvalue :',eigmna
        endif
      endif
      return
      end
