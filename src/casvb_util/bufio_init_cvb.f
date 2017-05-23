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
      subroutine bufio_init_cvb(file_id1)
      implicit real*8 (a-h,o-z)
      logical tstfile_cvb
      external tstfile_cvb
#include "bufio_cvb.fh"
#include "idbl_cvb.fh"

      file_id=file_id1
      ibuf=0
      if(.not.tstfile_cvb(file_id))then
        nbuf=0
        dnbuf=DBLE(nbuf)
        call wrlow_cvb(dnbuf,1,file_id,0)
      else
        call rdlow_cvb(dnbuf,1,file_id,0)
        nbuf=nint(dnbuf)
      endif
      nword=lbuf/idbl
      call izero(izbuffer,lbuf)
      return
      end
