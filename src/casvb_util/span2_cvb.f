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
      subroutine span2_cvb(c,nvec,s,n,metr)
      implicit real*8 (a-h,o-z)
#include "malloc_cvb.fh"
      common /span_comcvb/iaddr,nvecmx,nvtot
      dimension c(n,nvec),s(*)

      if(nvtot.ne.0)then
        call span_cvb(w(iaddr),nvtot,nvtot,s,n,metr)
        call fmove(w(iaddr),c,n*nvtot)
      endif
      nvec=nvtot
      call mfreer_cvb(iaddr)
      return
      end
