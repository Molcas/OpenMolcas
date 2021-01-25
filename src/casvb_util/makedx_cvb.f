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
      subroutine makedx_cvb(dx,nparm,ioptc,
     >  heigvec,heigval,dxp,gradp,w2,
     >  close2conv,converged,nposeig,scalesmall,
     >  wrongstat,nnegeig,opth,alfastart,alfa)
      implicit real*8 (a-h,o-z)
#include "formats_cvb.fh"
      logical opth,close2conv,converged,wrongstat,scalesmall
#include "locopt1_cvb.fh"
#include "locopt2_cvb.fh"
#include "tune_cvb.fh"
      dimension heigvec(nparm,nparm),heigval(nparm)
      dimension dxp(nparm),dx(nparm),gradp(nparm),w2(nparm)
      save zero,one
      data zero/0.d0/,one/1.d0/

      alfa=alfastart
c  << Update is small >>
      if(cnrm.lt.hh.and.scalesmall)then
c  Close to wrong stationary point :
        if(wrongstat)then
          if(dnrm2_(nparm,gradp,1).lt.grdwrngtol)then
            write(6,*)' Gradient too small - not using information!'
            call fzero(w2,nparm)
            do i=1,nnegeig
            if(heigval(i).ge.eigwrngtol)w2(i)=sign(one,gradp(i))
            enddo
            do i=nnegeig+1,nparm
            if(heigval(i).le.-eigwrngtol)w2(i)=sign(one,gradp(i))
            enddo
            call getdxp_cvb(dxp,w2,heigval,nnegeig,nparm,alfa)
            cnrm=dnrm2_(nparm,dxp,1)
          endif
          call dscal_(nparm,hh/cnrm,dxp,1)
          cnrm=hh
          goto 600
        endif
        if(.not.opth.and.ip.ge.2)write(6,form2AF)
     >     ' Scaling update from :',cnrm,' to :',hh
        call dscal_(nparm,hh/cnrm,dxp,1)
        cnrm=hh
        goto 600
      elseif(cnrm.lt.hh)then
        goto 600
      endif

      call optalf_cvb(heigval,gradp,nparm,hh,alfa,
     >  nnegeig,alfastart,alftol)
      call getdxp_cvb(dxp,gradp,heigval,nnegeig,nparm,alfa)
      call expec_cvb(dxp,gradp,heigval,nnegeig,nparm,exp,exp1,exp2)
      cnrm=dnrm2_(nparm,dxp,1)
      if(.not.opth.and.ip.ge.2)write(6,formAF)
     >   ' Alpha and norm of update :',alfa,cnrm

600   continue


      if(ioptc.gt.0.and.(.not.opth.and.cnrm.lt.cnrmtol))then
        if(ip.ge.0)then
          write(6,'(a)')' '
          write(6,formAD)' WARNING - predicted update too small :',
     >      cnrm,cnrmtol
        endif
        ioptc=-2
        return
      endif
1100  call expec_cvb(dxp,gradp,heigval,nnegeig,nparm,exp,exp1,exp2)
      if(exp1.lt.-exp12tol.or.exp2.gt.exp12tol)then
        call dscal_(nparm,0.9d0,dxp,1)
        cnrm=dnrm2_(nparm,dxp,1)
        if(cnrm.lt.cnrmtol)then
          write(6,formAD)' Norm of update too small :',cnrm,cnrmtol
          call abend_cvb()
        endif
        goto 1100
      endif
      if(ip.ge.2.and.close2conv.and.
     >  (exp1.lt.zero.or.exp2.gt.zero))then
        write(6,*)' Warning - not a max/min direction !'
        if(nnegeig.gt.0)write(6,*)' Expected change for maximized',
     >    ' variables :',exp1
        if(nposeig.gt.0)write(6,*)' Expected change for minimized',
     >    ' variables :',exp2
      endif
      call mxatb_cvb(heigvec,dxp,nparm,nparm,1,dx)
      return
c Avoid unused argument warnings
      if (.false.) call Unused_logical(converged)
      end
