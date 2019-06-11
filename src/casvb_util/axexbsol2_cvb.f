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
      subroutine axexbsol2_cvb(ap,rhsp,itdav,maxdav,nfrdim1,
     >  solp,solp_res,eig,eig_res,
     >  eigval,eigvec,dxp,gradp,w2)
c  Solve linear equation in Davidson subspace.
      implicit real*8 (a-h,o-z)
#include "direct_cvb.fh"
#include "locopt1_cvb.fh"
#include "tune_cvb.fh"
      dimension ap(maxdav,maxdav),rhsp(maxdav)
      dimension solp(maxdav),solp_res(maxdav)
      dimension eigval(itdav),eigvec(itdav,itdav)
      dimension dxp(itdav),gradp(itdav),w2(itdav)
      save zero,one
      data zero/0d0/,one/1d0/

      do 100 it=1,itdav
100   call fmove_cvb(ap(1,it),eigvec(1,it),itdav)

      if(ip.ge.3)then
        write(6,*)' AP matrix :'
        do i=1,itdav
        eigval(i)=eigvec(i,i)
        eigvec(i,i)=eigvec(i,i)+corenrg
        enddo
        call mxprintd_cvb(eigvec,itdav,itdav,0)
        do i=1,itdav
        eigvec(i,i)=eigval(i)
        enddo
        write(6,*)' RHSP vector :'
        call mxprint_cvb(rhsp,1,itdav,0)
      endif

      call mxdiag_cvb(eigvec,eigval,itdav)

      if(ip.ge.2)then
        write(6,'(a)')' Eigenvalues :'
        do i=1,itdav
        eigval(i)=eigval(i)+corenrg
        enddo
        call vecprint_cvb(eigval,itdav)
        do i=1,itdav
        eigval(i)=eigval(i)-corenrg
        enddo
      endif

      call mxatb_cvb(rhsp,eigvec,1,itdav,itdav,gradp)
      if(ifollow.eq.1)then
        nposeig=nroot-1
        nnegeig=itdav-nposeig
      elseif(ifollow.eq.2)then
        nnegeig=nroot-1
        nposeig=itdav-nnegeig
      else
        write(6,*)' Error in IFOLLOW with direct Fletcher!',ifollow
        call abend_cvb()
      endif

      eigmx=-one
      eigmn=one
      if(nnegeig.gt.0)eigmx=eigval(nnegeig)
      if(nposeig.gt.0)eigmn=eigval(nnegeig+1)
      safety_use=safety
200   continue
      if(eigmx.lt.-signtol.and.eigmn.gt.signtol)then
        alfastart=zero
      else
        alfastart=max(eigmx,-eigmn,zero)+safety_use
      endif
      call getdxp_cvb(dxp,gradp,eigval,nnegeig,itdav,alfastart)
      cnrm=dnrm2_(itdav,dxp,1)
      if(alfastart.ne.zero)then
        gnrm=dnrm2_(itdav,gradp,1)
        if(cnrm.gt.1d-15.and.gnrm.gt.1d-15.and.safety_use.ne.1d-4)then
          ovr_dx_grad=ddot_(itdav,dxp,1,gradp,1)/(cnrm*gnrm)
          if(ovr_dx_grad.lt..3d0)then
            safety_use=1d-4
            goto 200
          endif
        endif
      endif

      call makedx_cvb(solp,itdav,0,
     >  eigvec,eigval,dxp,gradp,w2,
     >  .false.,.false.,nposeig,.false.,
     >  .false.,nnegeig,.false.,alfastart,eig)

      eig_res=eig
      call fmove_cvb(solp,solp_res,itdav)
      if(ip.ge.2)then
        write(6,'(a,f15.8)')' Eigenvalue :',eig
        write(6,'(a)')' Solution vector :'
        call vecprint_cvb(solp,itdav)
      endif
      return
c Avoid unused argument warnings
      if (.false.) call Unused_integer(nfrdim1)
      end
