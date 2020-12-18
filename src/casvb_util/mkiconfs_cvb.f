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
      subroutine mkiconfs_cvb()
      implicit real*8 (a-h,o-z)
      logical need_cas
#include "ext_cvb.fh"
#include "main_cvb.fh"
#include "optze_cvb.fh"
#include "files_cvb.fh"
#include "print_cvb.fh"

#include "frag_cvb.fh"
#include "formats_cvb.fh"
#include "malloc_cvb.fh"

c  ICONFS
      call rdioff_cvb(4,recinp,ioffs)
      call rdis_cvb(iw(ll(15)),nconf*noe,recinp,ioffs)
      return
      entry mksymelm_cvb()
      call rdioff_cvb(8,recinp,ioffs)
      call rdr_cvb(w(ls(1)),nsyme*norb*norb,recinp,ioffs)
      if(ip(2).ge.1.and..not.up2date_cvb('PRSYMELM'))then
        do 300 isyme=1,nsyme
        write(6,'(/,a,i4,3x,a)')' Symmetry element no.',
     >    isyme,tags(isyme)
        ishft=norb*norb*(isyme-1)
        call mxprint_cvb(w(ishft+ls(1)),norb,norb,0)
300     continue
        if(nsyme.gt.0)write(6,*)' '
        call untouch_cvb('PRSYMELM')
      endif
      return
      entry mkconstruc_cvb()
      call construc_cvb(w(ls(15)),iw(ls(16)))
      return
      entry mkrdint_cvb()
      return
      entry mkrdcas_cvb()
      if(ifinish.eq.0)then
        need_cas=icrit.eq.1.or.projcas
      else
        need_cas=ifcasci_cvb().and.(.not.variat)
      endif
      if(.not.need_cas)return
c  Get CASSCF eigenvector
      if(.not.ifcasci_cvb())then
        if(ip(1).ge.0.and.valid_cvb(strtci))
     >    call prtfid_cvb(' Warning: CI vector not found - no ',
     >    strtci)
        if(icrit.eq.1)then
          write(6,*)' No optimization without CASSCF vector!'
          call abend_cvb()
        endif
      else
        if(ip(3).ge.2)
     >    write(6,'(/,a)') ' Read CASSCF eigenvector:'
        call getci_cvb(w(lc(1)))
      endif
      call cinorm2_cvb(w(lc(1)),cnrm)
      cnrm=one/cnrm
      call ciscale2_cvb(w(lc(1)),cnrm,iscf,cscf)
      if((.not.up2date_cvb('CASCHECK')).or.ip(3).ge.2)then
        call untouch_cvb('CASCHECK')
c  Some checks
        if(abs(cnrm-one).gt.1.d-3)then
          if(ip(3).ge.0)write(6,formE)
     >      ' WARNING: Norm of CI vector read differs from one :',cnrm
        elseif(ip(3).ge.2)then
          write(6,formE)' Norm of CI vector read ',cnrm
        endif
        if(ip(3).ge.2.and.iscf.ne.0)then
          write(6,'(a,i6)')' SCF determinant:',iscf
          write(6,formE) '     coefficient:',cscf
        endif
        if(ifhamil_cvb())then
          call cicopy_cvb(w(lc(1)),w(lc(2)))
          call applyh_cvb(w(lc(2)))
          call cidot_cvb(w(lc(1)),w(lc(2)),eexp)
          if(ip(3).ge.1)write(6,formE)' CASSCF energy :',eexp+corenrg
          if(ip(3).ge.1)write(6,'(a)')' '
        endif
      endif
      if(.not.memplenty)call ciwr_cvb(w(lc(1)),61001.2d0)
      return
      end
