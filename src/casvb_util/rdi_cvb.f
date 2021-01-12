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
      subroutine rdi_cvb(ivec,n,file_id,ioffset)
      implicit real*8 (a-h,o-z)
#include "idbl_cvb.fh"
      dimension ivec(n),ibuf(8)

      call rdi_cvb_internal(ivec,ibuf)
*
*     This is to allow type punning without an explicit interface
      contains
      subroutine rdi_cvb_internal(ivec,ibuf)
      use iso_c_binding
      integer, target :: ivec(*),ibuf(*)
      real*8, pointer :: vec(:),buf(:)
      nreals=n/idbl
      nrem=n-nreals*idbl
      call c_f_pointer(c_loc(ivec(1)),vec,[nreals])
      call rdlow_cvb(vec,nreals,file_id,ioffset)
      nullify(buf)
      if(nrem.gt.0)then
        call c_f_pointer(c_loc(ibuf(1)),buf,[1])
        call rdlow_cvb(buf,1,file_id,nreals+ioffset)
        nullify(buf)
        call imove_cvb(ibuf,ivec(1+nreals*idbl),nrem)
      endif
      return
      end subroutine rdi_cvb_internal
*
      end
