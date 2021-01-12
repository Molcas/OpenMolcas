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
      subroutine decl_cvb(chr)
      implicit real*8 (a-h,o-z)
      character*(*) chr
#include "make_cvb.fh"

      iobj=0
      do 100 i=1,nobj
      if(charobj(i).eq.chr)iobj=i
100   continue
      if(iobj.gt.0)then
        if(iprint.gt.1)write(6,*)' Make object exists already :',chr
        return
      endif
      nobj=nobj+1
      if(nobj.gt.mxobj)then
        write(6,*)' Too many make objects, max :',mxobj
        call abend_cvb()
      endif
      charobj(nobj)=chr
      up2date(nobj)=.false.
      ioffs(nobj+1)=ioffs(nobj)
      joffs(nobj+1)=joffs(nobj)
      if(iprint.ge.10)then
        write(6,*)' IOFFS :',(ioffs(ii),ii=1,nobj+1)
        write(6,*)' JOFFS :',(joffs(ii),ii=1,nobj+1)
      endif
      return
      end
