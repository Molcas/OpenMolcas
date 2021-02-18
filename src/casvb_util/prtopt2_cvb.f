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
      subroutine prtopt2_cvb(iopt1,ioptim,italter,noptim,
     >  iorts,ifxorb,ifxstr,idelstr)
      implicit real*8 (a-h,o-z)
#include "main_cvb.fh"
#include "optze_cvb.fh"
#include "files_cvb.fh"
#include "print_cvb.fh"

#include "spinb_cvb.fh"
#include "malloc_cvb.fh"
      parameter(nmeth=12)
      character*3 ayn
      character*8 methkw(nmeth)
      character*9 sbformat
      dimension iorts(2*norb*(norb-1)/2)
      dimension ifxorb(norb),ifxstr(nvb),idelstr(nvb)
      save methkw
      data methkw/'Fletcher','    TRIM','Trustopt','Davidson',
     >            '   Steep','  Vb2cas',' AugHess','AugHess2',
     >            '   Check',' dFletch','    None','Super-CI'/

      if(ifinish.eq.0)then
        if(ip(3).ge.1.or.(ip(3).eq.0.and.
     >    (iopt1.eq.0.or.italter.eq.1)))then
          if(noptim.eq.1)then
            write(6,'(/,a)')
     >        ' -- Starting optimization ------------------'
          else
            write(6,'(/,a,i3,a)')' -- Starting optimization - step',
     >      ioptim,' --------'
          endif
        endif
        if(ip(3).ge.1.and.(iopt1.eq.0.or.italter.eq.1))then
          if(icrit.eq.1)then
            write(6,'(/,a)')' Overlap-based optimization (Svb).'
          elseif(icrit.eq.2)then
            write(6,'(/,a)')' Energy-based optimization (Evb).'
          endif
          write(6,'(/,a,11x,a)')' Optimization algorithm:',
     >      methkw(imethod)
          write(6,'(a,i13)')' Maximum number of iterations:',mxiter
          ayn=' No'
          if(projcas)ayn='Yes'
          if(projcas)write(6,'(a,31x,a)')' Casproj:',ayn
          ayn=' No'
          if(projsym)ayn='Yes'
          if(projsym)write(6,'(a,31x,a)') ' Symproj:',ayn
          ayn=' No'
          sbformat='(a,19x,a)'
          write(sbformat(4:5),'(i2)')31-len_trim_cvb(spinb(kbasis))
          write(6,sbformat)' Spin basis:',
     >      spinb(kbasis)(1:len_trim_cvb(spinb(kbasis)))
          if(isaddle.gt.0)write(6,'(/,a,i9)')
     >      ' Saddle-point optimization, order:',isaddle
          if(nort.gt.0)then
            write(6,'(/,i4,a)')nort,' orthogonalization pairs defined :'
            write(6,6100)(ior,iorts(1+(ior-1)*2),
     >                        iorts(2+(ior-1)*2),ior=1,nort)
6100        format(3(i4,': ',i2,' -',i2))
          endif
          if(nfxorb.eq.norb)then
            write(6,'(/,a)')' All orbitals will be frozen.'
          elseif(nfxorb.gt.0)then
            write(6,'(/,a)')' Following orbitals will be frozen :'
            itmp = mstacki_cvb(nfxorb)
            ifx=0
            do 700 i=1,norb
            if(ifxorb(i).ge.0.and.ifxorb(i).le.norb)then
              ifx=ifx+1
              iw(ifx+itmp-1)=i
            endif
700         continue
            nfxorb=ifx
            write(6,'(14i3)')(iw(ii+itmp-1),ii=1,nfxorb)
            call mfreei_cvb(itmp)
          endif
          if(nfxvb.gt.0.and.lfxvb.eq.0)then
            write(6,'(/,a)')' Following structures will be frozen :'
            write(6,'(14i3)')(ifxstr(ii),ii=1,nfxvb)
          elseif(nfxvb.eq.0.and.lfxvb.eq.1)then
            write(6,'(/,a)')' All structures will be frozen.'
          elseif(nfxvb.gt.0.and.lfxvb.eq.1)then
            write(6,'(/,2a)')' Following structures coefficients',
     >        ' will be optimized :'
            write(6,'(14i3)')(ifxstr(ii),ii=1,nfxvb)
          endif
          if(nzrvb.gt.0.and.lzrvb.eq.0)then
            write(6,'(/,a)')' Following structures will be deleted :'
            write(6,'(14i3)')(idelstr(ii),ii=1,nzrvb)
          elseif(nzrvb.eq.0.and.lzrvb.eq.1)then
            write(6,'(/,a)')' All structures will be deleted.'
          elseif(nzrvb.gt.0.and.lzrvb.eq.1)then
            write(6,'(/,a)')' Following structures will not be',
     >        ' deleted :'
            write(6,'(14i3)')(idelstr(ii),ii=1,nzrvb)
          endif
          write(6,'(/,a)')
     >      ' -------------------------------------------'
        endif
        if(iopt1.eq.0.or.italter.eq.1)call tuneprint_cvb()
      elseif(ifinish.lt.3)then
        if(ip(3).ge.1.or.(ip(3).eq.0.and.
     >    (iopt1.eq.0.or.italter.eq.1)))then
          if(noptim.eq.1)then
            write(6,'(/,a)')
     >        ' -- Wavefunction summary -------------------'
          else
            write(6,'(/,a,i3,a)')' -- Wavefunction summary - step',
     >        ioptim,' ---------'
          endif
        endif
      endif
      return
      end
c  ******************************************
c  ** Highest-level input parsing routines **
c  ******************************************
