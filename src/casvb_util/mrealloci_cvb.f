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
      subroutine mrealloci_cvb(ipoint,nword)
      implicit real*8 (a-h,o-z)
#include "idbl_cvb.fh"
#include "memman_cvb.fh"

      iraddr=(ipoint-1)/idbl+1
      nwordr=(nword+idbl-1)/idbl
      call mreallocr_cvb(iraddr,nwordr)
      ipoint=(iraddr-1)*idbl+1
      if(memdebug)write(6,*)'   mrealloci : nword & pointer :',
     >  nword,ipoint
      return
      end
c
c  -- Zeroing routines - just front-ends ---
c
      integer function mheaprz_cvb(nword)
      implicit real*8 (a-h,o-z)
#include "memman_cvb.fh"
#include "malloc_cvb.fh"

      if(memdebug)write(6,*)' mheaprz :'
      mheaprz_cvb=mheapr_cvb(nword)
      call fzero(w(mheaprz_cvb),nword)
      return
      end
      integer function mstackrz_cvb(nword)
      implicit real*8 (a-h,o-z)
#include "memman_cvb.fh"
#include "malloc_cvb.fh"

      if(memdebug)write(6,*)' mstackrz :'
      mstackrz_cvb=mstackr_cvb(nword)
      call fzero(w(mstackrz_cvb),nword)
      return
      end
      integer function mheapiz_cvb(nword)
      implicit real*8 (a-h,o-z)
#include "memman_cvb.fh"
#include "malloc_cvb.fh"

      if(memdebug)write(6,*)' mheapiz :'
      mheapiz_cvb=mheapi_cvb(nword)
      call izero(iw(mheapiz_cvb),nword)
      return
      end
      integer function mstackiz_cvb(nword)
      implicit real*8 (a-h,o-z)
#include "memman_cvb.fh"
#include "malloc_cvb.fh"

      if(memdebug)write(6,*)' mstackiz :'
      mstackiz_cvb=mstacki_cvb(nword)
      call izero(iw(mstackiz_cvb),nword)
      return
      end
c  ****************************
c  ** Time and date routines **
c  ****************************
      function tim_cvb(cpu0)
      implicit real*8(a-h,o-z)
      call timing(cpu,cpusince,wall,wallsince)
      tim_cvb=cpu-cpu0
      return
      end
      function tim0_cvb()
      implicit real*8(a-h,o-z)
      call timing(cpu,cpusince,wall,wallsince)
      tim0_cvb=cpu
      return
      end
