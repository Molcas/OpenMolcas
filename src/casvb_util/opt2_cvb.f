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
      subroutine opt2_cvb(orbs,cvb)
      implicit real*8 (a-h,o-z)
#include "ext_cvb.fh"
#include "main_cvb.fh"
#include "optze_cvb.fh"
#include "files_cvb.fh"
#include "print_cvb.fh"

#include "formats_cvb.fh"
#include "malloc_cvb.fh"
      dimension orbs(norb,norb),cvb(nvb)

c The following initialization is to appease a compiler
      fx=zero
      ioptc=0
      iter=0

      if(imethod.eq.11)then
c  Method = None :
        goto 10
      endif
      if(imethod.eq.4)then
        if(icrit.eq.1)then
          call svbd_cvb(orbs,cvb,fx,ioptc,iter)
        elseif(icrit.eq.2)then
          call evbd_cvb(orbs,cvb,fx,ioptc,iter)
        endif
        goto 10
      elseif(imethod.eq.6)then
        call evb2cas_cvb(orbs,cvb,fx,ioptc,iter)
        goto 10
      endif

      fx=zero

      call optize_cvb(fx,ioptc,iter,
     >  imethod,isaddle,mxiter,icrit.eq.1,corenrg,ip(3),ip(4)-2,ip(4)-2,
     >  strucopt)

      if(ioptc.eq.-1.and.mxiter.gt.0)then
        if(ip(3).ge.0)write(6,'(a,i4)')
     >    ' Maximum number of iterations reached:',mxiter
        if(ip(3).ge.0)
     >  write(6,'(a)')' Calculation NOT converged!!!'
      endif
10    continue
      if(icrit.eq.1)then
        svb=fx
      else
        evb=fx
      endif
      if(ip(5).ge.0)then
        if(imethod.ne.11)then
          if(icrit.eq.1)write(6,formE)' Final Svb :',svb
          if(icrit.eq.2)write(6,formE)' Final Evb :',evb
        endif
        if(ip(3).le.1.and.ioptc.ne.-1)
     >    write(6,'(a,i4)')' Number of iterations used:',iter
      endif
      if(ip(5).ge.2)then
        call report_cvb(orbs,norb)
        write(6,'(/,a)')' Structure coefficients :'
        write(6,'(a)')' ------------------------'
        call vecprint_cvb(cvb,nvb)
      endif
      convinone=((ioptc.eq.0.and.iter.le.1).or.endvar)
      ioptc_new=ioptc
      n_iter=n_iter+iter
      if(ioptc.eq.1)ioptc_new=mxiter
      if(ioptc.eq.0)ioptc_new=iter
      return
      end
