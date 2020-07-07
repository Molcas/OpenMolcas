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
      subroutine bufio_chbuf_cvb(jbuf)
c  Change buffer position to JBUF. NB: Does not write current buffer,
c  neither is JBUF buffer read.
      implicit real*8 (a-h,o-z)
#include "bufio_cvb.fh"
#include "idbl_cvb.fh"

c  Dumy writes so that we don't exceed end-of-file:
      do 100 kbuf=nbuf+1,jbuf-1
      ibuf=kbuf
      call bufio_wrzbuf_cvb()
100   continue
      ibuf=jbuf
      return
      end
