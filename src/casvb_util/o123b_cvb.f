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
      subroutine o123b_cvb(nparm,
     >  dxnrm,grdnrm,close2conv)
      implicit real*8 (a-h,o-z)
      logical close2conv
#include "malloc_cvb.fh"
#include "opt2_cvb.fh"

      call o123b2_cvb(nparm,
     >  w(ix(1)),w(ix(3)),w(ix(4)),w(ix(5)),w(ix(6)),w(ix(7)),
     >  dxnrm)
      return
c Avoid unused argument warnings
      if (.false.) then
        call Unused_real(grdnrm)
        call Unused_logical(close2conv)
      end if
      end
