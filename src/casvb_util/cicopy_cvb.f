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
      subroutine cicopy_cvb(cvec1,cvec2)
      implicit real*8(a-h,o-z)
#include "ext_cvb.fh"
#include "main_cvb.fh"
#include "optze_cvb.fh"
#include "files_cvb.fh"
#include "print_cvb.fh"

#include "malloc_cvb.fh"
      dimension cvec1(*),cvec2(*)
c  *********************************************************************
c  *                                                                   *
c  *  CICOPY  := Copy CI vector                                        *
c  *                                                                   *
c  *********************************************************************

      ivec1=nint(cvec1(1))
      ivec2=nint(cvec2(1))
      iformat=iform_ci(ivec1)
      iform_ci(ivec2)=iform_ci(ivec1)
      call setcnt2_cvb(ivec2,igetcnt2_cvb(ivec1))
      if(iformat.eq.0)then
        call fmove(w(iaddr_ci(ivec1)),w(iaddr_ci(ivec2)),ndet)
      else
        write(6,*)' Unsupported format in CICOPY :',iformat
        call abend_cvb()
      endif
      return
      end
