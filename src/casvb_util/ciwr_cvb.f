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
      subroutine ciwr_cvb(cvec,recn)
      implicit real*8(a-h,o-z)
#include "main_cvb.fh"
#include "optze_cvb.fh"
#include "files_cvb.fh"
#include "print_cvb.fh"

#include "malloc_cvb.fh"
      dimension cvec(*)
c  *********************************************************************
c  *                                                                   *
c  *  CIWR   := Write CI vector.                                       *
c  *                                                                   *
c  *********************************************************************

      ivec=nint(cvec(1))
      iformat=iform_ci(ivec)
      if(iformat.eq.0)then
        ioffs=0
        call wris_cvb(iform_ci(ivec),1,recn,ioffs)
        call wris_cvb(icnt_ci(ivec),1,recn,ioffs)
        call wrrs_cvb(w(iaddr_ci(ivec)),ndet,recn,ioffs)
      else
        write(6,*)' Unsupported format in CIWR :',iformat
        call abend_cvb()
      endif
      return
      end
