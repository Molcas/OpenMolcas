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
      subroutine mkstrtgs_cvb(orbsao,irdorbs,cvb,recn,kbasis1)
      implicit real*8 (a-h,o-z)
      logical use_ao,ifmos_cvb
#include "ext_cvb.fh"
#include "main_cvb.fh"
#include "optze_cvb.fh"
#include "files_cvb.fh"
#include "print_cvb.fh"

#include "malloc_cvb.fh"
#include "mo_cvb.fh"
      dimension orbsao(nbas_mo,norb),irdorbs(norb),cvb(*)

      call rdheader_cvb(recn,norb1,nbas_mo1,nvb1,kbasiscvb,
     >  ioffs_orbs,ioffs_cvb,ioffs_orbsao,ioffs_orbslao)
      use_ao=(ifmos_cvb().and.((.not.variat).or.(variat.and.
     >  nmcscf.eq.1)).and.nbas_mo.eq.nbas_mo1.and.ioffs_orbsao.gt.0)
      do 100 iorb=1,norb
      if(.not.use_ao)then
        irdorbs(iorb)=1
        call rdgspr_cvb(recn,orbsao(1,iorb),iorb,norb,1,ierr)
      else
        irdorbs(iorb)=2
        call rdgspr_cvb(recn,orbsao(1,iorb),iorb,nbas_mo,3,ierr)
      endif
      if(ierr.ne.0)then
        call prtfid_cvb(' Error in orbital read from ',recn)
        write(6,'(a)')     ' Orbital no :',iorb
        if(use_ao)then
          write(6,'(a)')     ' AO basis ? : Yes'
        else
          write(6,'(a)')     ' AO basis ? : No'
        endif
        call abend_cvb()
      endif
100   continue
      call rdgspr_cvb(recn,cvb,1,nvb,2,ierr)
      return
c Avoid unused argument warnings
      if (.false.) call Unused_integer(kbasis1)
      end
