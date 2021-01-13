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
      subroutine mkafter_cvb(chr1,chr2)
      implicit real*8 (a-h,o-z)
      character*(*) chr1,chr2
#include "make_cvb.fh"

      call undepend2_cvb(chr1,chr2,1)

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
      ndep_ij=ndep_ij+1
      if(ndep_ij.gt.mxdep)then
        write(6,*)' Too many make dependencies, max :',mxdep
        call abend_cvb()
      endif
      do 200 i=ioffs(nobj+1),ioffs(iobj+1)+1,-1
      i_dep_on_j(i+1)=i_dep_on_j(i)
200   continue
      i_dep_on_j(ioffs(iobj+1)+1)=jobj
      do 300 i=iobj+1,nobj+1
      ioffs(i)=ioffs(i)+1
300   continue
      return
      end
