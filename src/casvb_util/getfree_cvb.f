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
      subroutine getfree_cvb(nfrr,n_div,nfrdim,iter,fx)
      implicit real*8 (a-h,o-z)
      logical orb_is_cheap
#include "ext_cvb.fh"
#include "main_cvb.fh"
#include "optze_cvb.fh"
#include "files_cvb.fh"
#include "print_cvb.fh"

#include "formats_cvb.fh"
#include "malloc_cvb.fh"
      save fxlast

      dxmove=.true.
      if(iter.ge.0)then
        if(ip(3).ge.2)then
          write(6,'(/,a,i5,a,f10.3,a)')' Iteration',iter,' at',
     >     tim_cvb(cpu0),' CPU seconds'
          write(6,'(a)')' ---------------------------------------'
        endif
        if(icrit.eq.1)then
          if(ip(3).ge.2)write(6,formE)' Svb :      ',fx
          if(ip(3).ge.2.and.iter.gt.1)
     >      write(6,formE)' Svb chg. : ',fx-fxlast
        elseif(icrit.eq.2)then
          if(ip(3).ge.2)write(6,formE)' Evb :      ',fx
          if(ip(3).ge.2.and.iter.gt.1)
     >      write(6,formE)' Evb chg. : ',fx-fxlast
        endif
        if(ip(3).ge.2)then
          call report_cvb(w(lv(1)),norb)
          if(strucopt)then
            write(6,'(/,a)')' Structure coefficients :'
            write(6,'(a)')' ------------------------'
            call vecprint_cvb(w(lv(2)),nvb)
          endif
        endif
      endif
      fxlast=fx
      call make_cvb('ORBFREE')
      call make_cvb('CIFREE')
      nfrr=nfr
      if(imethod.ne.4)then
        nfrdim=max(0,nfr-1)
      else
        nfrdim=nfr
      endif
      orb_is_cheap=(icrit.eq.1.and..not.(proj.or.projcas))
c  Set N_DIV :
      if((.not.strucopt).or.(.not.orb_is_cheap))then
        n_div=0
      else
        n_div=nfrorb
      endif
      return
      end
