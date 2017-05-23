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
      subroutine evb2cas2_cvb(orbs,cvb,ioptc,iter,fx,
     >   dxnrm,dx_amx,
     >   civec,civb,civbh,res,resh,
     >   cvbdet,gjorb)
      implicit real*8 (a-h,o-z)
#include "ext_cvb.fh"
#include "main_cvb.fh"
#include "optze_cvb.fh"
#include "files_cvb.fh"
#include "print_cvb.fh"

#include "formats_cvb.fh"
#include "tols_cvb.fh"
      logical dx_ok,grad_ok
      dimension orbs(norb,norb),cvb(nvb)
      dimension civec(ndet),civb(ndet),civbh(ndet)
      dimension cvbdet(ndetvb)
      dimension gjorb(*)
      dimension h(2,2),eig(2)

      if(ip(3).ge.0)then
        write(6,'(/,a)')' Starting VB2CAS optimization.'
        write(6,'(a)')' -----------------------------'
      endif

      dx_ok=(dx_amx.lt.dx(1,3).and.dxnrm.lt.dx(2,3))

      if(.not.projcas)then
        call str2vbc_cvb(cvb,cvbdet)
        call vb2cic_cvb(cvbdet,civb)
      elseif(projcas)then
        if(memplenty)then
          call cicopy_cvb(civec,civbh)
        else
          call cird_cvb(civbh,61001.2d0)
        endif
        call fmove(orbs,orbinv,norb*norb)
        call mxinv_cvb(orbinv,norb)
        call gaussj_cvb(orbinv,gjorb)
        call applyt_cvb(civbh,gjorb)
        call pvbcopy_cvb(civbh,civb)
        call ci2vbc_cvb(civbh,cvbdet)
      endif

      call gaussj_cvb(orbs,gjorb)
      call applyt_cvb(civb,gjorb)
      call proj_cvb(civb)

      call cinorm_cvb(civb,cnrm)
      call ciscale_cvb(civb,one/sqrt(cnrm))

      call cicopy_cvb(civb,civbh)
      call applyh_cvb(civbh)

      call cidot_cvb(civb,civbh,evb)
      if(ip(3).ge.2)write(6,formAF)
     >  ' Residual calculation based on Evb :',evb+corenrg
c RES()=CIVBH()-EVB*CIVB()
      call cicopy_cvb(civbh,res)
      call cidaxpy_cvb(-evb,civb,res)

      if(tstfile_cvb(67123.2d0))then
        call cird_cvb(resh,67123.2d0)
        call cidot_cvb(res,resh,rescas_ovr)
c        call cidot_cvb(civb,resh,civb_ovr)
c        dxnrm_ci=sqrt(2d0*(one-civb_ovr))
c        write(6,*)' dxnrm dxnrm_ci :',dxnrm,dxnrm_ci
c        write(6,*)' gradient in VB basis :',2d0*rescas_ovr/dxnrm
        grad_ok=(2d0*rescas_ovr/dxnrm.lt.grd(1,3))
      else
        grad_ok=.false.
      endif
      call ciwr_cvb(civb,67123.2d0)

      call cinorm_cvb(res,resnrm)
      if(ip(3).ge.2)then
        write(6,'(a)')' '
        write(6,formAD)' Residual norm:',resnrm
        write(6,'(a)')' '
      endif
      call ciscale_cvb(res,one/sqrt(resnrm))
      call cidot_cvb(res,civb,ovr)
c RES()=RES()-OVR*CIVB()
      call cidaxpy_cvb(-ovr,civb,res)
      call cinorm_cvb(res,resnrm)
      call ciscale_cvb(res,one/sqrt(resnrm))
      call cidot_cvb(civbh,civb,h(1,1))
      call cidot_cvb(civbh,res,h(1,2))

      call cicopy_cvb(res,resh)
      call applyh_cvb(resh)

      call cidot_cvb(resh,civb,h(2,1))
      call cidot_cvb(resh,res,h(2,2))

      if(ip(3).ge.2)then
        write(6,*)' 2x2 Hamiltonian matrix :'
        eig(1)=h(1,1)
        eig(2)=h(2,2)
        h(1,1)=h(1,1)+corenrg
        h(2,2)=h(2,2)+corenrg
        call mxprintd_cvb(h,2,2,0)
        h(1,1)=eig(1)
        h(2,2)=eig(2)
      endif

      call mxdiag_cvb(h,eig,2)
      if(ip(3).ge.2)then
        write(6,*)' Eigenvalues :',eig(1)+corenrg,eig(2)+corenrg
        write(6,*)' Eigenvectors :'
        call mxprint_cvb(h,2,2,0)
      endif

      if(abs(h(1,1)).gt.abs(h(1,2)))then
        if(ip(3).ge.2)write(6,*)' Using root 1 :'
        call ciscale_cvb(civb,h(1,1))
        call cidaxpy_cvb(h(2,1),res,civb)
      else
        if(ip(3).ge.2)write(6,*)' Using root 2 :'
        call ciscale_cvb(civb,h(1,2))
        call cidaxpy_cvb(h(2,2),res,civb)
      endif
      call cinorm_cvb(civb,cnrm)
      call ciscale_cvb(civb,one/sqrt(cnrm))
      if(memplenty)then
c        call cidot_cvb(civb,civec,ovr)
        call cicopy_cvb(civb,civec)
      else
        call cird_cvb(res,61001.2d0)
c        call cidot_cvb(civb,res,ovr)
        call ciwr_cvb(civb,61001.2d0)
      endif

      evb=evb+corenrg
      fx=evb
      ovraa=one
      iter=1
      ioptc=0
c      ovrcrit=.125d-10
c      if(abs(one-abs(ovr)).gt.ovrcrit)iter=2
      if(.not.(dx_ok.and.grad_ok))iter=2
      call setcnt_cvb(civec,1)
      return
      end
