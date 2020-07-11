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
      subroutine touchdepend_cvb(chr1,chr2)
      implicit real*8 (a-h,o-z)
      character*(*) chr1,chr2
#include "make_cvb.fh"

      call undepend2_cvb(chr1,chr2,2)

      iobj=0
      jobj=0
      do 100 i=1,nobj
      if(charobj(i).eq.chr1)iobj=i
      if(charobj(i).eq.chr2)jobj=i
100   continue
      if(iobj.eq.0)then
        write(6,*)' Make object not found :',chr1
        call abend_cvb()
      endif
      if(jobj.eq.0)then
        write(6,*)' Make object not found :',chr2
        call abend_cvb()
      endif
      ndep_ji=ndep_ji+1
      if(ndep_ji.gt.mxdep)then
        write(6,*)' Too many make dependencies, max :',mxdep
        call abend_cvb()
      endif
      do 200 i=joffs(nobj+1),joffs(jobj+1)+1,-1
      j_dep_on_i(i+1)=j_dep_on_i(i)
200   continue
      j_dep_on_i(joffs(jobj+1)+1)=iobj
      do 300 i=jobj+1,nobj+1
      joffs(i)=joffs(i)+1
300   continue

      if(.not.up2date(jobj))up2date(iobj)=.false.
      return
      end
