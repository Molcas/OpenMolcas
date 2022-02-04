*s***********************************************************************
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
      subroutine daname_cvb(lu,fname)
      implicit real*8 (a-h,o-z)
      character*(*) fname
#include "dancom_cvb.fh"
      integer find_lu, isfreeunit, i_open
      logical find_unused, is_opened
      external find_lu, isfreeunit, is_opened

      ilu=find_lu(fname)
      if(ilu.gt.0)then
        lu=ilu
        goto 1000
      endif
      find_unused=(lu.lt.1)
      if(.not.find_unused)find_unused=is_opened(lu)
      if(find_unused)then
        ilu=isfreeunit(10)
      endif

1000  if(is_opened(lu)) then
        i_open=1
      else
        i_open=0
        lu=10  ! initialize
      end if
      call istkpush_cvb(idan,i_open)
      if(i_open.eq.0)call daname(lu,fname)
      return
      end
      subroutine daclos_cvb(lu)
      implicit real*8 (a-h,o-z)
#include "dancom_cvb.fh"
c      logical find_unused

      call istkpop_cvb(idan,iwasopen)
      if(iwasopen.eq.0)call daclos(lu)
      return
      end
      subroutine daninit_cvb()
      implicit real*8 (a-h,o-z)
#include "dancom_cvb.fh"
c      logical find_unused
      call istkinit_cvb(idan,mxfiles)
      return
      end
