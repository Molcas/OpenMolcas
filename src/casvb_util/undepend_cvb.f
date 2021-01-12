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
      subroutine undepend_cvb(chr1,chr2)
      implicit real*8 (a-h,o-z)
      character*(*) chr1,chr2
#include "make_cvb.fh"

      ic=3
50    continue

      iobj=0
      jobj=0
      do 100 i=1,nobj
      if(charobj(i).eq.chr1)iobj=i
      if(charobj(i).eq.chr2)jobj=i
100   continue
      if(iobj.eq.0)then
        if(mustdeclare)then
          write(6,*)' Make object not found :',chr1
          call abend_cvb()
        endif
        call decl_cvb(chr1)
        goto 50
      endif
      if(jobj.eq.0)then
        if(mustdeclare)then
          write(6,*)' Make object not found :',chr2
          call abend_cvb()
        endif
        call decl_cvb(chr2)
        goto 50
      endif

      if(iprint.ge.10)write(6,*)' Cancel I depends on J :',iobj,jobj
      n_cancelled=0
      if(mod(ic,2).eq.1)then
190     continue
        do 200 i=ioffs(iobj)+1,ioffs(iobj+1)
        if(i_dep_on_j(i).eq.jobj)then
          do 300 j=i,ioffs(nobj+1)-1
          i_dep_on_j(j)=i_dep_on_j(j+1)
300       continue
          do 400 ii=iobj+1,nobj+1
          ioffs(ii)=ioffs(ii)-1
400       continue
          n_cancelled=n_cancelled+1
          goto 190
        endif
200     continue
      endif

      m_cancelled=0
      if(ic.ge.2)then
490     continue
        do 500 i=joffs(jobj)+1,joffs(jobj+1)
        if(j_dep_on_i(i).eq.iobj)then
          do 600 j=i,joffs(nobj+1)-1
          j_dep_on_i(j)=j_dep_on_i(j+1)
600       continue
          do 700 ii=jobj+1,nobj+1
          joffs(ii)=joffs(ii)-1
700       continue
          m_cancelled=m_cancelled+1
          goto 490
        endif
500     continue
      endif

      ndep_ij=ndep_ij-n_cancelled
      ndep_ji=ndep_ji-m_cancelled
      return
      end

      subroutine undepend2_cvb(chr1,chr2,ic1)
      implicit real*8 (a-h,o-z)
      character*(*) chr1,chr2
#include "make_cvb.fh"
      ic=ic1
50    continue

      iobj=0
      jobj=0
      do 100 i=1,nobj
      if(charobj(i).eq.chr1)iobj=i
      if(charobj(i).eq.chr2)jobj=i
100   continue
      if(iobj.eq.0)then
        if(mustdeclare)then
          write(6,*)' Make object not found :',chr1
          call abend_cvb()
        endif
        call decl_cvb(chr1)
        goto 50
      endif
      if(jobj.eq.0)then
        if(mustdeclare)then
          write(6,*)' Make object not found :',chr2
          call abend_cvb()
        endif
        call decl_cvb(chr2)
        goto 50
      endif

      if(iprint.ge.10)write(6,*)' Cancel I depends on J :',iobj,jobj
      n_cancelled=0
      if(mod(ic,2).eq.1)then
190     continue
        do 200 i=ioffs(iobj)+1,ioffs(iobj+1)
        if(i_dep_on_j(i).eq.jobj)then
          do 300 j=i,ioffs(nobj+1)-1
          i_dep_on_j(j)=i_dep_on_j(j+1)
300       continue
          do 400 ii=iobj+1,nobj+1
          ioffs(ii)=ioffs(ii)-1
400       continue
          n_cancelled=n_cancelled+1
          goto 190
        endif
200     continue
      endif

      m_cancelled=0
      if(ic.ge.2)then
490     continue
        do 500 i=joffs(jobj)+1,joffs(jobj+1)
        if(j_dep_on_i(i).eq.iobj)then
          do 600 j=i,joffs(nobj+1)-1
          j_dep_on_i(j)=j_dep_on_i(j+1)
600       continue
          do 700 ii=jobj+1,nobj+1
          joffs(ii)=joffs(ii)-1
700       continue
          m_cancelled=m_cancelled+1
          goto 490
        endif
500     continue
      endif

      ndep_ij=ndep_ij-n_cancelled
      ndep_ji=ndep_ji-m_cancelled
      return
      end
