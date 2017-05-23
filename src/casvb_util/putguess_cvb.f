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
      subroutine putguess_cvb(orbs,cvb,recn)
      implicit real*8 (a-h,o-z)
      logical use_ao,ifmos_cvb
#include "ext_cvb.fh"
#include "main_cvb.fh"
#include "optze_cvb.fh"
#include "files_cvb.fh"
#include "print_cvb.fh"

#include "malloc_cvb.fh"
#include "mo_cvb.fh"
      dimension orbs(norb,*),cvb(*)

      call wrheader_cvb(recn,norb,nbas_mo,nvb,kbasiscvb,
     >  ioffs_orbs,ioffs_cvb,ioffs_orbsao,ioffs_orbslao)
      call rdheader_cvb(recn,norb1,nbas_mo1,nvb1,kbasiscvb1,
     >  ioffs_orbs,ioffs_cvb,ioffs_orbsao,ioffs_orbslao)
      do 100 iorb=1,norb
100   call wrgspr_cvb(recn,orbs(1,iorb),iorb,norb,1,ierr)
      call wrgspr_cvb(recn,cvb,1,nvb,2,ierr)
      use_ao=(ifmos_cvb().and.((.not.variat).or.(variat.and.
     >  endvar)))
      if(use_ao)then
        iorbsao = mstackr_cvb(nbas_mo*norb)
        call mo2ao_cvb(orbs,w(iorbsao),norb)
        do 200 iorb=1,norb
200     call wrgspr_cvb(recn,w((iorb-1)*nbas_mo+iorbsao),
     >    iorb,nbas_mo,3,ierr)
        if(ip(5).ge.2)then
          write(6,'(/,a)')' VB orbitals in AO basis :'
          write(6,'(a)')  ' -------------------------'
          call mxprint_cvb(w(iorbsao),nbas_mo,norb,0)
        endif
        if(ploc)then
          i1 = mstackr_cvb(norb*norb)
          i2 = mstackr_cvb(norb*norb)
          i3 = mstackr_cvb(norb)
          call getr_plc(w(i1))
          call transp_cvb(w(i1),w(i1),norb,norb)
          call mxatb_cvb(w(i1),orbs,norb,norb,norb,w(i2))
          call lmo2ao_cvb(w(i2),w(iorbsao),norb)
          do 300 iorb=1,norb
300       call wrgspr_cvb(recn,w((iorb-1)*nbas_mo+iorbsao),
     >      iorb,nbas_mo,4,ierr)
          if(ip(5).ge.2)then
            write(6,'(/,a)')
     >        ' Original localized VB orbitals in AO basis :'
            write(6,'(a)')
     >        ' --------------------------------------------'
            call mxprint_cvb(w(iorbsao),nbas_mo,norb,0)
          endif
          do 400 i=1,norb
          w(i+i3-1)=dnrm2_(norb,w((i-1)*norb+i2),1)
400       call dscal_(norb,1d0/w(i+i3-1),w((i-1)*norb+i2),1)
          if(ip(5).ge.2)then
            write(6,'(/,a)')' Norms of original localized VB orbitals :'
            write(6,'(a)')  ' -----------------------------------------'
            call mxprint_cvb(w(i3),1,norb,0)
          endif
          call mfreer_cvb(i1)
        endif
        call mfreer_cvb(iorbsao)
      endif
      return
      end
