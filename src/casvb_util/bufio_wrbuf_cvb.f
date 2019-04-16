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

      call bufio_wrbuf_internal(ibuffer)
      return

      entry bufio_wrzbuf_cvb()
      call bufio_wrzbuf_internal(izbuffer)
      return

      entry bufio_rdbuf_cvb()
      call bufio_rdbuf_internal(ibuffer)
      return
*
*     This is to allow type punning without an explicit interface
      contains
      subroutine bufio_wrbuf_internal(ibuffer)
      use iso_c_binding
      integer, target :: ibuffer(*)
      real*8, pointer :: buffer(:)
      if(ibuf.eq.0)return

      ioffset=(ibuf-1)*lbuf/idbl
      call c_f_pointer(c_loc(ibuffer(1)),buffer,[nword])
      call wrlow_cvb(buffer,nword,file_id,ioffset+1)
      nullify(buffer)
      if(ibuf.gt.nbuf)then
        nbuf=ibuf
      endif
      return
      end subroutine bufio_wrbuf_internal
      subroutine bufio_wrzbuf_internal(izbuffer)
      use iso_c_binding
      integer, target :: izbuffer(*)
      real*8, pointer :: buffer(:)
      if(ibuf.eq.0)return

      ioffset=(ibuf-1)*lbuf/idbl
      call c_f_pointer(c_loc(izbuffer(1)),buffer,[nword])
      call wrlow_cvb(buffer,nword,file_id,ioffset+1)
      nullify(buffer)
      if(ibuf.gt.nbuf)then
        nbuf=ibuf
      endif
      return
      end subroutine bufio_wrzbuf_internal
      subroutine bufio_rdbuf_internal(ibuffer)
      use iso_c_binding
      integer, target :: ibuffer(*)
      real*8, pointer :: buffer(:)
      if(nbuf.lt.ibuf)then
        call izero(ibuffer,lbuf)
        return
      endif
      ioffset=(ibuf-1)*lbuf/idbl
      call c_f_pointer(c_loc(ibuffer(1)),buffer,[nword])
      call rdlow_cvb(buffer,nword,file_id,ioffset+1)
      nullify(buffer)
      return
      end subroutine bufio_rdbuf_internal
*
      end
