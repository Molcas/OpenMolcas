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
      subroutine ciscale2_cvb(cvec,scale,iscf,cscf)
      implicit real*8(a-h,o-z)
#include "ext_cvb.fh"
#include "main_cvb.fh"
#include "optze_cvb.fh"
#include "files_cvb.fh"
#include "print_cvb.fh"

#include "malloc_cvb.fh"
      dimension cvec(*)

      ivec=nint(cvec(1))
      iscf=0
      cscf=zero
      iformat=iform_ci(ivec)
      if(iformat.eq.0)then
        do 100 idet=1,ndet
        w(idet+iaddr_ci(ivec)-1)=scale*w(idet+iaddr_ci(ivec)-1)
        if(abs(w(idet+iaddr_ci(ivec)-1)).gt.p8)then
          iscf=idet
          cscf=w(idet+iaddr_ci(ivec)-1)
        endif
100     continue
      else
        write(6,*)' Unsupported format in CISCALE2 :',iformat
        call abend_cvb()
      endif
      return
      end
