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
      subroutine trust_cvb(iopth,opth,maxize,fx,fxbest,exp,
     >  hh,dxnrm,ioptc,scalesmall1,close2conv,converged,skipupd)
      implicit real*8 (a-h,o-z)
#include "formats_cvb.fh"
      logical opth,maxize,scalesmall1,close2conv,converged
      logical skipupd
      logical dfx_ok,zz_ok
#include "trst_cvb.fh"
#include "tune_cvb.fh"
#include "print_cvb.fh"
      save zero,half,one
      data zero/0d0/,half/.5d0/,one/1d0/

      call zz_cvb(dum,zz,fx,fxbest,exp,-1)
      skipupd=.false.
      if(.not.close2conv)then
        ipu=1
      else
        ipu=2
      endif
      nopth=nopth1(ipu)+nopth2(ipu)
      scalesmall1=scalesmall(ipu)

c
c---- Trust region control, based on the quadratic model
c
      iop=mod(iopth,nopth)
      if(iop.eq.0)iop=iop+nopth
      if(iopth.gt.0)cpropt(iop)=fx
c  Set HHOPT to HH actually used (might have been scaled ) :
      if(iopth.gt.0)hhopt(iop)=dxnrm
1     continue
      if(iopth.gt.0.and.mod(iopth,nopth).eq.0.and.opth)then
c  Optimisation of trust region completed
        opth=.false.
        if(maxize)then
          call findmx_cvb(cpropt,nopth,cprbst,icprbst)
        else
          call findmn_cvb(cpropt,nopth,cprbst,icprbst)
        endif
        dfx=cprbst-fxbest
        dfx_ok=((dfx.gt.dfxmin(ipu).and.maxize).or.
     >          (dfx.lt.-dfxmin(ipu).and..not.maxize))
        zz_ok=(zz.gt.zzrejmin(ipu).and.zz.lt.zzrejmax(ipu))
        if(dfx_ok.and.zz_ok)then
c  << Accepting update >>
          iopth=0
c  Restore HH as used before (so repeated update vector will be
c  exactly the same) :
          if(icprbst.le.nopth1(ipu))then
            hh=hhkeep*(one+(DBLE(icprbst)-half*DBLE(nopth1(ipu)+1))
     >        *delopth1(ipu))
          elseif(icprbst.le.nopth)then
            ! IFG: nopth1 was used in these two calls, probably a bug
            if(maxize)then
              call findmx_cvb(cpropt,nopth,cprbst,icprbst2)
            else
              call findmn_cvb(cpropt,nopth,cprbst,icprbst2)
            endif
            hh=hhkeep*(one+(DBLE(icprbst2)-half*DBLE(nopth1(ipu)+1))
     >        *delopth1(ipu))
            gap2=hhkeep*delopth1(ipu)*delopth2(ipu)
            hh=hh+gap2*
     >        (DBLE(icprbst-nopth1(ipu))-half*DBLE(nopth2(ipu)+1))
          endif
          hh=min(hh,hhmax(ipu))
c  Scale trust region size according to ZZ :
          if(zz.lt.zzacclim(1,ipu))then
            scale=hhaccfac(1,ipu)
          elseif(zz.gt.zzacclim(4,ipu))then
            scale=hhaccfac(5,ipu)
          elseif(zz.lt.zzacclim(2,ipu))then
            scale=hhaccfac(2,ipu)
          elseif(zz.gt.zzacclim(3,ipu))then
            scale=hhaccfac(4,ipu)
          else
            scale=hhaccfac(3,ipu)
          endif
          if(nopth.gt.1)then
            if(hhopt(icprbst).gt.1d-8.and.hh/hhopt(icprbst).gt.2d0)then
              hhkeep=hh
            else
c  Scale trust region size according to ZZ :
              hhkeep=hh*scale
            endif
          else
            oldstep=hhopt(icprbst)
c  Scale trust region size according to ZZ :
            if(scale.ge.one)then
              hhkeep=max(hhkeep,oldstep*scale)
            else
              hhkeep=hhkeep*scale
            endif
          endif
          skipupd=(icprbst.eq.nopth)
        else
c  << Rejecting update >>
          if(converged)then
            iopth=0
            hh=zero
            return
          endif
          if(ip(3).ge.1)write(6,'(a)')' Rejecting step.'
          call findmn_cvb(hhopt,nopth,hh_min,idum)
          hhkeep=min(hh_min,hhkeep)*hhrejfac(ipu)
          gap2=hhkeep*delopth1(ipu)*delopth2(ipu)
          if(nopth2(ipu).eq.0)gap2=zero
          hhlargest=hhkeep*(one+(DBLE(nopth1(ipu))
     >      -half*DBLE(nopth1(ipu)+1))*delopth1(ipu))
     >      +gap2*(DBLE(nopth-nopth1(ipu))-half*DBLE(nopth2(ipu)+1))
          if(hhlargest.lt.hhtol(ipu))then
            if(ip(3).ge.0)then
              write(6,formAD)
     >          ' Trust region size smaller than tolerance !',
     >          hhlargest,hhtol(ipu)
              write(6,'(a)')' Calculation NOT converged!'
            endif
            ioptc=-2
            return
          endif
          goto 1
        endif
      else
        if(iopth.eq.0.and.nopth.gt.1.and.(ip(3).ge.2))
     >    write(6,'(/,a)')' Optimising trust region size :'
        opth=.true.
        iopth=iopth+1
        ioptst=mod(iopth,nopth)
        if(ioptst.eq.0)ioptst=ioptst+nopth

        if(ioptst.le.nopth1(ipu))then
          hh=hhkeep*(one+(DBLE(ioptst)-half*DBLE(nopth1(ipu)+1))
     >      *delopth1(ipu))
        elseif(ioptst.le.nopth)then
          ! IFG: nopth1 was used in these two calls, probably a bug
          if(maxize)then
            call findmx_cvb(cpropt,nopth,cprbst,icprbst)
          else
            call findmn_cvb(cpropt,nopth,cprbst,icprbst)
          endif
          hh=hhkeep*(one+(DBLE(icprbst)-half*DBLE(nopth1(ipu)+1))
     >      *delopth1(ipu))
          gap2=hhkeep*delopth1(ipu)*delopth2(ipu)
          hh=hh+gap2*(DBLE(ioptst-nopth1(ipu))
     >         -half*DBLE(nopth2(ipu)+1))
        endif
        hh=min(hh,hhmax(ipu))
        hhopt(ioptst)=hh
      endif
      return
      end
