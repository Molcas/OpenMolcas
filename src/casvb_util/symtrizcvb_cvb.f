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
      subroutine symtrizcvb_cvb(vecstr)
      implicit real*8 (a-h,o-z)
#include "main_cvb.fh"
#include "optze_cvb.fh"
#include "files_cvb.fh"
#include "print_cvb.fh"

#include "WrkSpc.fh"
      dimension vecstr(nvb)
      dimension dum(1)

      if(iconstruc.eq.0)then
        return
      elseif(iconstruc.eq.1)then
        i1 = mstackr_cvb(ndetvb)
        i2 = mstackr_cvb(nvb)
        call symtrizcvb2_cvb(vecstr,
     >    iwork(ls(13)),iwork(ls(16)),work(i1),work(i2))
        call mfreer_cvb(i1)
        call symtrizcvb3_cvb(vecstr,iwork(ls(10)))
      elseif(iconstruc.eq.2)then
        call schmidtd_cvb(work(ls(15)),nconstr,vecstr,1,dum,nvb,0)
      endif
      return
      end
