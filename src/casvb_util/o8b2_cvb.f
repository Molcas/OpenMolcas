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
      subroutine o8b2_cvb(nparm,
     >  dx,grad,eigvec,eigval,
     >  dxnrm,close2conv)
      implicit real*8 (a-h,o-z)
      logical close2conv
#include "opt_cvb.fh"
#include "locopt1_cvb.fh"
#include "locopt2_cvb.fh"
#include "trst_cvb.fh"
#include "tune_cvb.fh"

      dimension dx(nparm),grad(nparm)
      dimension eigvec(nparm+1,nparm+1),eigval(nparm+1)
      save one
      data one/1d0/

      call fzero(eigvec,(nparm+1)*(nparm+1))
      do 100 iprm=1,nparm
      eigvec(iprm+1,1)=grad(iprm)
      eigvec(1,iprm+1)=grad(iprm)
      eigvec(iprm+1,iprm+1)=one
      call hess_cvb(eigvec(2,iprm+1))
100   continue
      write(6,*)' Augmented Hessian matrix :'
      call mxprint_cvb(eigvec,nparm+1,nparm+1,0)
      call mxdiag_cvb(eigvec,eigval,nparm+1)
      iroot=nparm+1
      if(ip.ge.2)then
        write(6,'(a)')' Eigenvalues of augmented Hessian :'
        call vecprint_cvb(eigval,nparm+1)
        write(6,'(a)')' Eigenvector to be followed :'
        call vecprint_cvb(eigvec(1,iroot),nparm+1)
      endif
      write(6,*)' Following root no :',iroot
      call fmove_cvb(eigvec(2,iroot),dx,nparm)
      if(abs(eigvec(1,iroot)).gt.1d-8)then
        fac1=one/eigvec(1,iroot)
      else
        fac1=sign(one,eigvec(1,iroot))
      endif
      call dscal_(nparm,fac1,dx,1)
      dxnrm=dnrm2_(nparm,dx,1)
      if(.not.close2conv)then
        ipu=1
      else
        ipu=2
      endif
      if(dxnrm.gt.hh.or.scalesmall(ipu))then
        call dscal_(nparm,hh/dxnrm,dx,1)
        dxnrm=hh
      endif
      return
      end
