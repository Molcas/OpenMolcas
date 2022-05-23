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
      subroutine untouch_cvb(chr)
      implicit real*8 (a-h,o-z)
      character*(*) chr
#include "make_cvb.fh"

50    iobj=0
      do 100 i=1,nobj
      if(charobj(i).eq.chr)iobj=i
100   continue
      if(iobj.eq.0)then
        if(mustdeclare)then
          write(6,*)' Make object not found :',chr
          call abend_cvb()
        endif
        call decl_cvb(chr)
        goto 50
      endif
      if(up2date(iobj))return
      if(iprint.ge.1)
     >  write(6,'(/,a,i3,2a)')' Untouch object no.',iobj,', name : ',
     >  charobj(iobj)
      up2date(iobj)=.true.
      return
      end
