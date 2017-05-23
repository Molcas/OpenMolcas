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
      subroutine prtopt_cvb()
      implicit real*8 (a-h,o-z)
#include "ext_cvb.fh"
#include "main_cvb.fh"
#include "optze_cvb.fh"
#include "files_cvb.fh"
#include "print_cvb.fh"

#include "seth_cvb.fh"
#include "initopt_cvb.fh"
#include "loopcntr_cvb.fh"
#include "malloc_cvb.fh"
      external istkprobe_cvb
      logical istkprobe_cvb

c  First determine if end of multi-step optimization may have been reached:
      if(istkprobe_cvb(istackrep))then
        call istkpop_cvb(istackrep,nc_zeroed)
        call istkpop_cvb(istackrep,nconvinone)
        call istkpop_cvb(istackrep,italter)
        call istkpop_cvb(istackrep,mxalter)
        call istkpop_cvb(istackrep,kk2)
        call istkpop_cvb(istackrep,ioptstep2)
        call istkpop_cvb(istackrep,ioptstep1)
        call istkpush_cvb(istackrep,ioptstep1)
        call istkpush_cvb(istackrep,ioptstep2)
        call istkpush_cvb(istackrep,kk2)
        call istkpush_cvb(istackrep,mxalter)
        call istkpush_cvb(istackrep,italter)
        call istkpush_cvb(istackrep,nconvinone)
        call istkpush_cvb(istackrep,nc_zeroed)
      else
        ioptstep1=0
        italter=0
      endif

      i1 = mstacki_cvb(2*norb*(norb-1)/2)
      i2 = mstacki_cvb(norb)
      i3 = mstacki_cvb(nvb)
      i4 = mstacki_cvb(nvb)

      call rdioff_cvb(11,recinp,ioffs)
      call rdis_cvb(iw(i2),norb,recinp,ioffs)
      call rdis_cvb(iw(i3),nfxvb,recinp,ioffs)
      call rdis_cvb(iw(i4),nzrvb,recinp,ioffs)
      call rdis_cvb(iw(i1),2*nort,recinp,ioffs)

      call prtopt2_cvb(ioptstep1,ioptim,italter,noptim,
     >  iw(i1),iw(i2),iw(i3),iw(i4))
      call mfreei_cvb(i1)
      return
      end
