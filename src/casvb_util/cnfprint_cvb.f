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
      subroutine cnfprint_cvb()
      implicit real*8 (a-h,o-z)
c ... Make: up to date? ...
      logical, external :: up2date_cvb
#include "main_cvb.fh"
#include "optze_cvb.fh"
#include "files_cvb.fh"
#include "print_cvb.fh"

#include "inpmod_cvb.fh"
#include "frag_cvb.fh"
#include "malloc_cvb.fh"
      logical recinpcmp_cvb
      dimension idum(1)

      if(recinpcmp_cvb(4))call touch_cvb('CNFPRINT')

      if(ip(1).ge.0.and.(.not.up2date_cvb('CNFPRINT')))then
        i1 = mstacki_cvb(max(noe,noe*nconf))
        call rdioff_cvb(1,recinp,ioffs)
        call rdis_cvb(idum,1,recinp,ioffs)
        !noe1=idum(1)
        call rdis_cvb(idum,1,recinp,ioffs)
        !nconf1=idum(1)
        call rdis_cvb(idum,1,recinp,ioffs)
        !kbasiscvb1=idum(1)
        call rdis_cvb(iw(i1),noe*nconf,recinp,ioffs)
        if(nconf.eq.0)then
          do 225 i=1,min(nel,norb)
          iw(i+i1-1)=1
225       continue
          do 250 i=1,nel-norb
          iw(i+i1-1)=2
250       continue
        endif
        nconf_off=0
        do 300 ifrag=1,nfrag
        if(nfrag.gt.1)write(6,'(/,a,i3)')
     >    ' Configuration list for wavefunction fragment',ifrag
        write(6,'(/,a)')' Spatial VB configurations'
        write(6,'(a)')' -------------------------'
        write(6,'(a)')'     Conf. =>   Orbitals'
        call cnfprt_cvb(iw(noe*nconf_off+i1),
     >    nconf_fr(ifrag),nel_fr(ifrag))
        write(6,'(/,a,i6)')' Number of VB configurations :',
     >    nconf_fr(ifrag)
        write(6,'(a,i6)')  '           VB structures     :',
     >    nvbr_fr(ifrag)
        write(6,'(a,i6)')  '           VB determinants   :',
     >    ndetvb_fr(ifrag)
        nconf_off=nconf_off+nconf_fr(ifrag)
300     continue
        call mfreei_cvb(i1)
        call make_cvb('CNFPRINT')
      endif
      return
      end
