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
      subroutine bufio_wrbuf_cvb()
      implicit real*8 (a-h,o-z)
#include "bufio_cvb.fh"
#include "idbl_cvb.fh"

      if(ibuf.eq.0)return

      ioffset=(ibuf-1)*lbuf/idbl
      call wrlow_cvb(ibuffer,nword,file_id,ioffset+1)
      if(ibuf.gt.nbuf)then
        nbuf=ibuf
      endif
      return

      entry bufio_wrzbuf_cvb()
      if(ibuf.eq.0)return

      ioffset=(ibuf-1)*lbuf/idbl
      call wrlow_cvb(izbuffer,nword,file_id,ioffset+1)
      if(ibuf.gt.nbuf)then
        nbuf=ibuf
      endif
      return

      entry bufio_rdbuf_cvb()
      if(nbuf.lt.ibuf)then
        call izero(ibuffer,lbuf)
        return
      endif
      ioffset=(ibuf-1)*lbuf/idbl
      call rdlow_cvb(ibuffer,nword,file_id,ioffset+1)
      return
      end
