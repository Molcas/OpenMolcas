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

      call wri_cvb_internal(ivec,ibuf)
*
*     This is to allow type punning without an explicit interface
      contains
      subroutine wri_cvb_internal(ivec,ibuf)
      use iso_c_binding
      integer, target :: ivec(*),ibuf(*)
      real*8, pointer :: vec(:),buf(:)
      nreals=n/idbl
      nrem=n-nreals*idbl
      if(nreals.gt.0) then
        call c_f_pointer(c_loc(ivec(1)),vec,[nreals])
        call wrlow_cvb(vec,nreals,file_id,ioffset)
        nullify(vec)
      end if
      if(nrem.gt.0)then
        len=0
        call lendat_cvb(file_id,len)
        if(len.ge.1+nreals+ioffset) then
          call c_f_pointer(c_loc(ibuf(1)),buf,[1])
          call rdlow_cvb(buf,1,file_id,nreals+ioffset)
          nullify(buf)
        end if
        call imove_cvb(ivec(1+nreals*idbl),ibuf,nrem)
c  Trying for a "clean" write (unwritten integers written as zeros),
c  but explicit zeroing of ibuf is not necessary as long as idbl <= 2.
c        call izero(ibuf(nrem+1),idbl-nrem)
        call c_f_pointer(c_loc(ibuf(1)),buf,[1])
        call wrlow_cvb(buf,1,file_id,nreals+ioffset)
        nullify(buf)
      endif
      return
      end subroutine wri_cvb_internal
*
      end
