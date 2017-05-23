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
      logical function up2date_cvb(chr)
      implicit real*8 (a-h,o-z)
      character*(*) chr
#include "make_cvb.fh"

      iobj=0
      do 100 i=1,nobj
100   if(charobj(i).eq.chr)iobj=i
      if(iobj.eq.0)then
        write(6,*)' Make object not found :',chr
        call abend_cvb()
      endif
      up2date_cvb=up2date(iobj)
      return
      end
