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
      subroutine o123b2_cvb(nparm,
     >  dx,eigvec,eigval,dxp,gradp,wrk,
     >  dxnrm)
      implicit real*8 (a-h,o-z)
#include "opt_cvb.fh"
#include "locopt1_cvb.fh"
#include "locopt2_cvb.fh"
#include "trst_cvb.fh"
#include "tune_cvb.fh"

      dimension dx(nparm)
      dimension eigvec(nparm,nparm),eigval(nparm)
      dimension dxp(nparm),gradp(nparm),wrk(nparm)
      save zero,one
      data zero/0.d0/,one/1d0/

      if(maxize)then
        nposeig=min(isaddle,nparm)
      else
        nposeig=max(nparm-isaddle,nparm)
      endif
      nnegeig=nparm-nposeig

      eigmx=-one
      eigmn=one
      if(nnegeig.gt.0)eigmx=eigval(nnegeig)
      if(nposeig.gt.0)eigmn=eigval(nnegeig+1)
      safety_use=safety
c200   continue
      if(eigmx.lt.-signtol.and.eigmn.gt.signtol)then
        alfastart=zero
      else
        alfastart=max(eigmx,-eigmn,zero)+safety_use
      endif
      call getdxp_cvb(dxp,gradp,eigval,nnegeig,nparm,alfastart)
      cnrm=dnrm2_(nparm,dxp,1)
c  Increased level shift not necessary when full Hessian is calculated :
c      if(alfastart.ne.zero)then
c        gnrm=dnrm2_(nparm,gradp,1)
c        if(cnrm.gt.1d-15.and.gnrm.gt.1d-15.and.safety_use.ne.1d-4)then
c          ovr_dx_grad=ddot_(nparm,dxp,1,gradp,1)/(cnrm*gnrm)
c          if(ovr_dx_grad.lt..3d0)then
c            safety_use=1d-4
c            goto 200
c          endif
c        endif
c      endif

      call makedx_cvb(dx,nparm,0,
     >  eigvec,eigval,dxp,gradp,wrk,
     >  .false.,.false.,nposeig,.false.,
     >  .false.,nnegeig,.false.,alfastart,eig)
      dxnrm=dnrm2_(nparm,dx,1)

      return
      end
