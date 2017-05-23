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
      subroutine scorr2_cvb(cvbdet,dvbdet,evbdet,
     >      ssq,wvbdet,iperm)
      implicit real*8 (a-h,o-z)
#include "ext_cvb.fh"
#include "main_cvb.fh"
#include "optze_cvb.fh"
#include "files_cvb.fh"
#include "print_cvb.fh"


#include "formats_cvb.fh"
      dimension cvbdet(ndetvb),dvbdet(ndetvb),evbdet(ndetvb)
      dimension ssq(norb,norb),wvbdet(ndetvb),iperm(norb)
      save cut,half
      data cut/1.d-10/,half/.5d0/

      write(6,'(/,1x,a)') 'Expectation values of (s(i)+s(j))**2'
      snorm=ddot_(ndetvb,cvbdet,1,dvbdet,1)
      ssnorm=ddot_(ndetvb,cvbdet,1,evbdet,1)
      write(6,formAF)
     >  ' Lower triangle uses SPIN function with Snorm=',ssnorm
      write(6,formAF)
     >  ' Upper triangle uses FULL function with Snorm=',snorm
c! DLC
c     snorm=one/snorm
c     ssnorm=one/ssnorm
      phase=(-one)**iabs(nalf-nbet)
      snorm=phase/snorm
      ssnorm=phase/ssnorm
c! DLC
      call fzero(ssq,norb*norb)
      tot=zero
      stot=zero
      do 100 mu=1,norb
      do 100 nu=mu+1,norb
c  Apply s_mu x s_nu to the wavefunction
      do 200 i=1,norb
200   iperm(i)=i
      iperm(mu)=nu
      iperm(nu)=mu
      call fmove(cvbdet,wvbdet,ndetvb)
      call permvb_cvb(wvbdet,iperm)
      sum=one-ddot_(ndetvb,wvbdet,1,dvbdet,1)*snorm
      ssum=one-ddot_(ndetvb,wvbdet,1,evbdet,1)*ssnorm
      tot=tot+sum
      stot=stot+ssum
      ssq(mu,nu)=sum
100   ssq(nu,mu)=ssum
      call mxprint_cvb(ssq,norb,norb,0)
      tot=tot+r3by4*dble(norb-2*norb*(norb-1)/2)
      stot=stot+r3by4*dble(norb-2*norb*(norb-1)/2)
      scheck=half*dble(iabs(nalf-nbet))*(half*dble(iabs(nalf-nbet))+one)
      if(abs(tot-scheck).gt.cut .or. abs(stot-scheck).gt.cut)
     >  write(6,formAD) 'WARNING: spins ',stot,tot,scheck
      return
      end
