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
      subroutine getmo_cvb(cmo,ic,ic2)
      implicit real*8 (a-h,o-z)
#include "mo_cvb.fh"
#include "malloc_cvb.fh"
      dimension cmo(*)

      i1 = mstackr_cvb(nbasisq_mo)

      if(ic.le.1)then
        i2 = mstackr_cvb(0)
        call getmo2_cvb(cmo,w(i2),w(i1),ic,ic2)
      else
        i2 = mstackr_cvb(nbas_mo*nbas_mo)
        call getmo2_cvb(w(i2),cmo,w(i1),ic,ic2)
      endif
      call mfreer_cvb(i1)
      return
      end
