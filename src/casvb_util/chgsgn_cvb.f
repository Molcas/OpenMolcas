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
      subroutine chgsgn_cvb(fx)
      implicit real*8 (a-h,o-z)
#include "main_cvb.fh"
#include "optze_cvb.fh"
#include "files_cvb.fh"
#include "print_cvb.fh"

#include "frag_cvb.fh"
#include "malloc_cvb.fh"

      if(nfrag.le.1)then
        call dscal_(nvb,-one,w(lv(2)),1)
        call dscal_(ndetvb,-one,w(lv(5)),1)
      else
        call dscal_(nvb_fr(1),-one,w(lv(2)),1)
        call dscal_(ndetvb_fr(1),-one,w(lv(5)),1)
      endif
      call touch_cvb('CVB')
      call fx_cvb(fx,.false.)
      return
      end
