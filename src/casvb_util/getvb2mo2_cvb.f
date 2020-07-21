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
      subroutine getvb2mo2_cvb(orbs,norbcheck,nelcheck)
      implicit real*8 (a-h,o-z)
#include "files_cvb.fh"
#include "casvb.fh"
      dimension orbs(*)

      call cvbinit_cvb()
      call rdheader_cvb(recn_vbwfn,norb1,nbas_mo1,nvb1,kbasiscvb1,
     >  ioffs_orbs,ioffs_cvb,ioffs_orbsao,ioffs_orbslao)
      if(norbcheck.ne.norb1)then
        call prtfid_cvb(' Error - present number of orbitals'//
     >    ' not consistent with number on ',recn_vbwfn)
        write(6,*)' Numbers :',norbcheck,norb1
        call abend_cvb()
      endif
      ioff=1
      do 100 iorb=1,norb1
      call  rdgspr_cvb(recn_vbwfn,orbs(ioff),iorb,norb1,1,ierr)
      if(ierr.ne.0)then
        write(6,*)' Error in VB orbital read :',ierr
        call abend()
      endif
      ioff=ioff+norb1
100   continue
      return
c Avoid unused argument warnings
      if (.false.) call Unused_integer(nelcheck)
      end
