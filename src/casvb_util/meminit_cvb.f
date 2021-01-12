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
      subroutine meminit_cvb()
      implicit real*8 (a-h,o-z)
#include "memman_cvb.fh"

      memdebug=.false.
      nfield=0
      ioff_r=0
      ioff_i=0
      call setmem('trace=off')
      call setmem('clear=off')
      if(memdebug)then
        write(6,*)' Casvb memory handler initialized.'
        write(6,*)' Memory offsets : integer= ',ioff_i,' real= ',ioff_r
        write(6,*)' No. of fields in use :',nfield
      endif

      return
      end
c
c  -- All memory allocated via real*8 routines ---
c
      integer function mheapr_cvb(nword1)
c  Memory allocator (heap). Returns pointer to NWORD real*8 words.
      implicit real*8 (a-h,o-z)
#include "memman_cvb.fh"

      nword=nword1
      if(memdebug)write(6,*)'     Enter mheapr: nword :',nword
      if(nword.lt.0)then
        write(6,*)' Error: attempting to allocate negative ',
     >    'amount of memory.'
        write(6,*)' nword=',nword
        call abend_cvb()
      endif
      call getmem('casvb','ALLO','REAL',ipoint_g,nword)
      mheapr_cvb=ipoint_g+ioff_r
      if(memdebug)write(6,*)'     mheapr: nword & pointer :',
     >  nword,mheapr_cvb
      return
      end
      integer function mstackr_cvb(nword)
c  Memory allocator (stack). Returns pointer to NWORD real*8 words.
      implicit real*8 (a-h,o-z)
#include "memman_cvb.fh"

      if(memdebug)write(6,*)'     Enter mstackr: nword :',nword
      mstackr_cvb=mheapr_cvb(nword)
      nfield=nfield+1
      if(nfield.gt.mxfield)then
        write(6,*)' Too many field in mstackr :',nfield,mxfield
        call abend_cvb()
      endif
      iaddr(nfield)=mstackr_cvb
      if(memdebug)write(6,*)'     mstackr: nword & pointer :',
     >  nword,mstackr_cvb
      return
      end
