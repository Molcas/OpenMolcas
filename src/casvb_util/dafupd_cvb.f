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
      subroutine dafupd_cvb(lu,ioffset)
      implicit real*8(a-h,o-z)
#include "fio.fh"
#include "idbl_cvb.fh"
      dimension ibuf(1000)
      Integer   mxddr
      data ibuf/1000*0/

      mxddr=1000
      nwrite=1000

      Call iDaFile(lu,8,ibuf,nwrite,mxddr)

      if(mxddr.lt.ioffset)then
        ioff=mxddr
100     nwrite=min((ioffset-ioff)*idbl,1000)
        call iDaFile(lu,1,ibuf,nwrite,ioff)
        if(ioff.lt.ioffset)goto 100
      endif
      return
      end
