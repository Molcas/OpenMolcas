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
      subroutine ddrestv_cvb(vec,avec,svec,ndim,ioffs,ause,suse)
      implicit real*8(a-h,o-z)
      logical ause,suse
#include "malloc_cvb.fh"
#include "direct_cvb.fh"
      dimension vec(ndim),avec(ndim),svec(ndim)

      nvguess=nvguess+1
      nvrestart=nvrestart+1
      if(nvguess.gt.maxd.or.nvrestart.gt.maxd)then
        write(6,*)' Too many guess vectors in Davidson!',nvguess,
     >    nvrestart,maxd
        call abend_cvb()
      endif
      if(ndim+ioffs.gt.nparm)then
        write(6,*)' Illegal call to DDRESTV :',ndim,ioffs,nparm
        call abend_cvb()
      endif
      iddvec=1
      call fzero(w(idd(iddvec)+(nvrestart-1)*nparm),ioffs)
      call fmove_cvb(vec,w(ioffs+idd(iddvec)+(nvrestart-1)*nparm),ndim)
      call fzero(w(ndim+ioffs+idd(iddvec)+(nvrestart-1)*nparm),
     >  nparm-ioffs-ndim)
      if(ause)then
        iddvec=iddvec+1
        call fzero(w(idd(iddvec)+(nvrestart-1)*nparm),ioffs)
        call fmove_cvb(avec,w(ioffs+idd(iddvec)+(nvrestart-1)*nparm),
     >    ndim)
        call fzero(w(ndim+ioffs+idd(iddvec)+(nvrestart-1)*nparm),
     >    nparm-ioffs-ndim)
      endif
      if(suse)then
        iddvec=iddvec+1
        call fzero(w(idd(iddvec)+(nvrestart-1)*nparm),ioffs)
        call fmove_cvb(svec,w(ioffs+idd(iddvec)+(nvrestart-1)*nparm),
     >    ndim)
        call fzero(w(ndim+ioffs+idd(iddvec)+(nvrestart-1)*nparm),
     >    nparm-ioffs-ndim)
      endif
      return
      end
