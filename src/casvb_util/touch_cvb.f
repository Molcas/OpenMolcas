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
      subroutine touch_cvb(chr)
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
      up2date(iobj)=.false.
      if(iprint.ge.1)
     >  write(6,'(/,a,i3,2a)')' Touch (1) of object no.',iobj,
     >  ', name : ',charobj(iobj)

c  Mark all "child" objects as out-of-date :
200   n_touched=0
      do 300 iobj=1,nobj
      if(.not.up2date(iobj))then
        do 400 i=joffs(iobj)+1,joffs(iobj+1)
        call touchrules_cvb(charobj(j_dep_on_i(i)))
        if(up2date(j_dep_on_i(i)))then
          up2date(j_dep_on_i(i))=.false.
          if(iprint.ge.1)
     >      write(6,'(/,a,i3,2a)')' Touch (2) of object no.',
     >      j_dep_on_i(i),', name : ',charobj(j_dep_on_i(i))
          n_touched=n_touched+1
        endif
400     continue
      endif
300   continue
      if(n_touched.ne.0)goto 200

      return
      end
