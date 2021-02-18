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
      subroutine daname_cvb(lu,fname)
      implicit real*8 (a-h,o-z)
      character*(*) fname
      parameter (mxfiles=40)
      common /dancom_cvb/idan(mxfiles)
#include "fio.fh"
      logical find_unused

      do 100 ilu=1,99
      if(isOpen(ilu).eq.1)then
        if(LuName(ilu).eq.fname)then
          lu=ilu
          goto 1000
        endif
      endif
100   continue
      find_unused=(lu.lt.1.or.lu.gt.99)
      if(.not.find_unused)find_unused=(isOpen(lu).eq.1)
      if(find_unused)then
        do 200 ilu=10,99
        if(isOpen(ilu).eq.0)then
          lu=ilu
          goto 1000
        endif
200     continue
        do 300 ilu=1,9
        if(ilu.ne.5.and.ilu.ne.6.and.isOpen(ilu).eq.0)then
          lu=ilu
          goto 1000
        endif
300     continue
        write(6,'(a)')' Unused unit number not found in DANAME!'
        call abend_cvb()
      endif

1000  call istkpush_cvb(idan,isOpen(lu))
      if(isOpen(lu).eq.0)call daname(lu,fname)
      return
      end
      subroutine daclos_cvb(lu)
      implicit real*8 (a-h,o-z)
      parameter (mxfiles=40)
      common /dancom_cvb/idan(mxfiles)
#include "fio.fh"
c      logical find_unused

      call istkpop_cvb(idan,iwasopen)
      if(iwasopen.eq.0)call daclos(lu)
      return
      end
      subroutine daninit_cvb()
      implicit real*8 (a-h,o-z)
      parameter (mxfiles=40)
      common /dancom_cvb/idan(mxfiles)
#include "fio.fh"
c      logical find_unused
      call istkinit_cvb(idan,mxfiles)
      return
      end
