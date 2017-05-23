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
      subroutine asonc10_cvb(c,axc,dum1,nvec,nprm)
      implicit real*8 (a-h,o-z)
#include "ext_cvb.fh"
#include "main_cvb.fh"
#include "optze_cvb.fh"
#include "files_cvb.fh"
#include "print_cvb.fh"

#include "malloc_cvb.fh"
      common /ipp/ipp,iter
      dimension c(nprm,nvec),axc(nprm,nvec)
c      save iter,ipp

      iter=iter+1
      if(ipp.ge.2)then
        write(6,'(/,a,i5,a,f10.3,a)')' Davidson iteration',iter,
     >    ' at',tim_cvb(cpu0),' CPU seconds'
        write(6,'(a)')
     >    ' -----------------------------------------------'
      endif

      do 100 ivec=1,nvec
      call fmove(c(1,ivec),axc(1,ivec),nprm)
      call hess_cvb(axc(1,ivec))
      call ddproj_cvb(axc(1,ivec),nprm)
100   continue
      return
c Avoid unused argument warnings
      if (.false.) call Unused_real(dum1)
      end

      subroutine asonc10init_cvb(ippinp)
      implicit real*8 (a-h,o-z)
#include "ext_cvb.fh"
#include "main_cvb.fh"
#include "optze_cvb.fh"
#include "files_cvb.fh"
#include "print_cvb.fh"

#include "malloc_cvb.fh"
      common /ipp/ ipp,iter
c      save iter

      iter=0
      ipp=ippinp
      call orthcvb_init_cvb()
      return
      end
