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
      subroutine ciscale_cvb(cvec,scale)
      implicit real*8(a-h,o-z)
#include "ext_cvb.fh"
#include "main_cvb.fh"
#include "optze_cvb.fh"
#include "files_cvb.fh"
#include "print_cvb.fh"

#include "malloc_cvb.fh"
      dimension cvec(*)
c  *********************************************************************
c  *                                                                   *
c  *  CISCALE := Analogous to the blas routine DSCAL                   *
c  *                                                                   *
c  *********************************************************************

      ivec=nint(cvec(1))
      iformat=iform_ci(ivec)
      if(iformat.eq.0)then
        call dscal_(ndet,scale,w(iaddr_ci(ivec)),1)
      else
        write(6,*)' Unsupported format in CISCALE :',iformat
        call abend_cvb()
      endif
      return
      end
