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
      subroutine free2all_cvb(vecfrom,vecto,nvec)
      implicit real*8 (a-h,o-z)
#include "ext_cvb.fh"
#include "main_cvb.fh"
#include "optze_cvb.fh"
#include "files_cvb.fh"
#include "print_cvb.fh"

#include "malloc_cvb.fh"
      dimension vecfrom(nfr,nvec),vecto(npr,nvec)

      do 100 ivec=1,nvec
      if(.not.orbfr_is_unit)then
        call mxatb_cvb(w(ls(14)),vecfrom(1,ivec),
     >    nprorb,nfrorb,1,vecto(1,ivec))
      else
        if(nprorb.gt.0) call fmove_cvb(vecfrom(1,ivec),
     >    vecto(1,ivec),nprorb)
      endif
      if(nprvb.gt.0) call fmove_cvb(vecfrom(nfrorb+1,ivec),
     >  vecto(nprorb+1,ivec),nprvb)
100   continue
      return
      end
