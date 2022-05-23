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
      subroutine ddsol72_cvb(hp,eigval,eigvec,dum,itdav,maxdav,nfrdim1,
     >  solp,solp_res,eig,eig_res)
c  Solve linear equation in Davidson subspace.
      implicit real*8 (a-h,o-z)
#include "direct_cvb.fh"
      dimension hp(maxdav,maxdav),eigval(itdav),eigvec(itdav,itdav)
      dimension solp(maxdav),solp_res(maxdav)

      if(ip.ge.3)then
        write(6,*)' HP matrix (b) :'
        call mxprint2_cvb(hp,itdav,maxdav,itdav,0)
      endif

      do it=1,itdav
      call fmove_cvb(hp(1,it),eigvec(1,it),itdav)
      enddo
      call mxdiag_cvb(eigvec,eigval,itdav)

      if(ifollow.le.2)then
        iroot=nroot
        jroot=mod(itdav,nroot)
        if(jroot.eq.0)jroot=nroot
        if(itdav.eq.maxdav.or.itdav.eq.nfrdim)jroot=nroot
        iroot=min(itdav,iroot)
        jroot=min(itdav,jroot)
        if(ifollow.eq.1)then
          iroot=itdav-iroot+1
          jroot=itdav-jroot+1
        endif
      elseif(ifollow.eq.3)then
        write(6,*)' Overlap-based root following not yet implemented!'
        call abend_cvb()
      elseif(ifollow.eq.4)then
c  Eigenvalue-based root following -- determine closest root :
        iroot=1
        delmin=abs(eigval(1)-eig)
        do 100 i=1,min(itdav,nroot)
        del=abs(eigval(i)-eig)
        if(del.lt.delmin)then
          delmin=del
          iroot=i
        endif
100     continue
        jroot=iroot
      endif
      eig=eigval(iroot)
      call fmove_cvb(eigvec(1,iroot),solp,itdav)
      eig_res=eigval(jroot)
      call fmove_cvb(eigvec(1,jroot),solp_res,itdav)
      if(ip.ge.2)then
        write(6,'(a)')' Eigenvalues :'
        call vecprint_cvb(eigval,itdav)
        write(6,'(a,i3,a)')' Eigenvector number',iroot,' :'
        call vecprint_cvb(solp,itdav)
        if(jroot.ne.iroot)then
          write(6,'(a,i3,a)')' Eigenvector number',jroot,' :'
        call vecprint_cvb(solp_res,itdav)
        endif
      endif
      return
c Avoid unused argument warnings
      if (.false.) then
        call Unused_real(dum)
        call Unused_integer(nfrdim1)
      end if
      end
      subroutine ddinit7_cvb(ifollow1,isaddle1,ip1)
      implicit real*8 (a-h,o-z)
#include "direct_cvb.fh"

      ifollow=ifollow1
      isaddle=isaddle1
      nroot=max(1,isaddle+1)
      ip=ip1
      return
      end
