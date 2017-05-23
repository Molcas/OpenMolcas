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
      subroutine make_cvb(chr)
      implicit real*8 (a-h,o-z)
      character*(*) chr
#include "make_cvb.fh"

50    iobj=0
      do 100 i=1,nobj
100   if(charobj(i).eq.chr)iobj=i
      if(iobj.eq.0)then
        if(mustdeclare)then
          write(6,*)' Make object not found :',chr
          call abend_cvb()
        endif
        call decl_cvb(chr)
        goto 50
      endif

c  Make sure all "parent" objects are up-to-date :
200   mkobj=iobj
300   continue
      do 400 i=ioffs(mkobj)+1,ioffs(mkobj+1)
      if(.not.up2date(i_dep_on_j(i)))then
        mkobj=i_dep_on_j(i)
        goto 300
      endif
400   continue
      if(.not.up2date(mkobj))then
        if(iprint.ge.1)
     >    write(6,'(/,a,i3,2a)')' Making object no.',mkobj,', name : ',
     >    charobj(mkobj)
        call rules_cvb(charobj(mkobj))
        up2date(mkobj)=.true.
      endif
      if(mkobj.ne.iobj)goto 200
      return
      end
