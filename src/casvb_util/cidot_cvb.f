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
      subroutine cidot_cvb(cvec1,cvec2,ret)
      implicit real*8(a-h,o-z)
#include "main_cvb.fh"
#include "optze_cvb.fh"
#include "files_cvb.fh"
#include "print_cvb.fh"

#include "malloc_cvb.fh"
      dimension cvec1(*),cvec2(*)
c  *********************************************************************
c  *                                                                   *
c  *  CIDOT  := Analogous to the blas routine DDOT                     *
c  *                                                                   *
c  *********************************************************************

      ivec1=nint(cvec1(1))
      ivec2=nint(cvec2(1))
      iformat1=iform_ci(ivec1)
      iformat2=iform_ci(ivec2)
      if(iformat1.ne.iformat2)then
        write(6,*)' Format discrepancy in CIDOT :',iformat1,iformat2
        call abend_cvb()
      endif
      if(iformat1.eq.0)then
        ret=ddot_(ndet,w(iaddr_ci(ivec1)),1,w(iaddr_ci(ivec2)),1)
      else
        write(6,*)' Unsupported format in CIDOT :',iformat1
        call abend_cvb()
      endif
      return
      end
