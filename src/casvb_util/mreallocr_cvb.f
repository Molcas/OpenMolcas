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
      subroutine mreallocr_cvb(ipoint,nword)
c  Memory allocator (heap). Reallocate pointer.
      implicit real*8 (a-h,o-z)
#include "memman_cvb.fh"
#include "malloc_cvb.fh"
#include "files_cvb.fh"

      if(memdebug)write(6,*)'     Enter mreallocr: nword & pointer :',
     >  nword,ipoint

      ipoint_g=ipoint-ioff_r

c      call getmem('casvb','CHAN','REAL',ipoint_g,nword)

c  Read and write data -- not efficient but safe and simple :
      call getmem('casvb','LENG','REAL',ipoint_g,nword_old)
      nword_move=min(nword,nword_old)
      call wrr_cvb(w(ipoint),nword_move,recn_tmp06,0)
      call mfreer_cvb(ipoint)
      ipoint=mheapr_cvb(nword)
      call rdr_cvb(w(ipoint),nword_move,recn_tmp06,0)

      if(memdebug)write(6,*)'     mreallocr : nword & pointer :',
     >  nword,ipoint
      return
      end
      integer function mavailr_cvb()
c  Memory allocator (heap). Returns number of real*8 words free.
      implicit real*8 (a-h,o-z)
#include "memman_cvb.fh"

      call getmem('casvb','MAX ','REAL',ipoint_g,nword)
      mavailr_cvb=nword
      if(memdebug)write(6,*)'     mavailr :',mavailr_cvb
      return
      end
c
c  -- Integer routines - just front-ends for real*8 ---
c
      integer function mheapi_cvb(nword)
      implicit real*8 (a-h,o-z)
#include "idbl_cvb.fh"
#include "memman_cvb.fh"

      if(memdebug)write(6,*)'   Enter mheapi: nword :',nword
      nwordr=(nword+idbl-1)/idbl
      iraddr=mheapr_cvb(nwordr)
      mheapi_cvb=(iraddr-1)*idbl+1
      if(memdebug)write(6,*)'   mheapi: nword & pointer :',
     >  nword,mheapi_cvb
      return
      end
      integer function mstacki_cvb(nword)
      implicit real*8 (a-h,o-z)
#include "idbl_cvb.fh"
#include "memman_cvb.fh"

      if(memdebug)write(6,*)'   Enter mstacki: nword :',nword
      nwordr=(nword+idbl-1)/idbl
      iraddr=mstackr_cvb(nwordr)
      mstacki_cvb=(iraddr-1)*idbl+1
      if(memdebug)write(6,*)'   mstacki: nword & pointer :',
     >  nword,mstacki_cvb
      return
      end
      integer function mavaili_cvb()
      implicit real*8 (a-h,o-z)
#include "idbl_cvb.fh"
#include "memman_cvb.fh"

      mavaili_cvb=mavailr_cvb()*idbl
      if(memdebug)write(6,*)'   mavaili :',mavaili_cvb
      return
      end
