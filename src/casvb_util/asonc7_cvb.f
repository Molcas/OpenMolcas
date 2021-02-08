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
      subroutine asonc7_cvb(c,axc,dum1,nvec,nprm)
      implicit real*8 (a-h,o-z)
#include "main_cvb.fh"
#include "optze_cvb.fh"
#include "files_cvb.fh"
#include "print_cvb.fh"

#include "malloc_cvb.fh"
      dimension c(nprm,nvec),axc(nprm,nvec)
      save igrad,iter,ipp
      save thresh
      data thresh/1d-15/

      iter=iter+1
      if(ipp.ge.2)then
        write(6,'(/,a,i5,a,f10.3,a)')' Davidson iteration',iter,
     >    ' at',tim_cvb(cpu0),' CPU seconds'
        write(6,'(a)')
     >    ' -----------------------------------------------'
      endif

      do 100 ivec=1,nvec
      axc(1,ivec)=ddot_(nprm-1,w(igrad),1,c(2,ivec),1)
      call fmove_cvb(c(2,ivec),axc(2,ivec),nprm-1)
c  Save Hessian application (& DNRM2 call) whenever possible :
c  (C assumed to be normalized)
      if(abs(abs(c(1,ivec))-one).gt.thresh)then
        call hess_cvb(axc(2,ivec))
      elseif(dnrm2_(nprm-1,axc(2,ivec),1).gt.thresh)then
        call hess_cvb(axc(2,ivec))
      endif
      call daxpy_(nprm-1,c(1,ivec),w(igrad),1,axc(2,ivec),1)
      call ddproj_cvb(axc(2,ivec),nprm-1)
100   continue
      return
      entry asonc7init_cvb(igradinp,ippinp)
      iter=0
      igrad=igradinp
      ipp=ippinp
      call orthcvb_init_cvb()
      return
c Avoid unused argument warnings
      if (.false.) call Unused_real(dum1)
      end
