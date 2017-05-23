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
      subroutine wri_cvb(ivec,n,file_id,ioffset)
      implicit real*8 (a-h,o-z)
#include "idbl_cvb.fh"
      dimension ivec(n)
      save ibuf
      dimension ibuf(8)
      data ibuf/8*0/

      nreals=n/idbl
      nrem=n-nreals*idbl
      if(nreals.gt.0)call wrlow_cvb(ivec,nreals,file_id,ioffset)
      if(nrem.gt.0)then
        len=0
        call lendat_cvb(file_id,len)
        if(len.ge.1+nreals+ioffset)
     >    call rdlow_cvb(ibuf,1,file_id,nreals+ioffset)
        call imove_cvb(ivec(1+nreals*idbl),ibuf,nrem)
c  Trying for a "clean" write (unwritten integers written as zeros),
c  but explicit zeroing of ibuf is not necessary as long as idbl <= 2.
c        call izero(ibuf(nrem+1),idbl-nrem)
        call wrlow_cvb(ibuf,1,file_id,nreals+ioffset)
      endif
      return
      end
