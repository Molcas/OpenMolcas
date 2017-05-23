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
      subroutine mkgrd_cvb(civb,civb2,grad,dvbdet,np,doorb)
      implicit real*8 (a-h,o-z)
      logical doorb
#include "ext_cvb.fh"
#include "main_cvb.fh"
#include "optze_cvb.fh"
#include "files_cvb.fh"
#include "print_cvb.fh"

#include "malloc_cvb.fh"
      dimension civb(ndet),civb2(ndet),grad(npr),dvbdet(ndetvb)

      call fzero(grad,nprorb)
      if(doorb)call onedens_cvb(civb,civb2,grad,.false.,1)
      if(strucopt)then
        call ci2vbg_cvb(civb2,dvbdet)
        if(np-nprorb.eq.nvb)then
          call vb2strg_cvb(dvbdet,grad(nprorb+1))
        elseif(np-nprorb.lt.nvb)then
          i1 = mstackrz_cvb(nvb)
          call vb2strg_cvb(dvbdet,w(i1))
          call fmove(w(i1),w(lv(5)),np-nprorb)
          call mfreer_cvb(i1)
        else
          write(6,*)' Error in mkgrd - np-nprorb > nvb :',np,nprorb,nvb
        endif
      endif
      return
      end
