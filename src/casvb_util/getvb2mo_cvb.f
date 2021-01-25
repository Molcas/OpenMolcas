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
      subroutine getvb2mo_cvb(orbs)
      implicit real*8 (a-h,o-z)
#include "files_cvb.fh"
#include "casvb.fh"
      dimension orbs(*)

      if(ifvb.eq.1)call cvbinit_cvb()
      call rdheader_cvb(recn_vbwfn,norb1,nbas_mo1,nvb1,kbasiscvb1,
     >  ioffs_orbs,ioffs_cvb,ioffs_orbsao,ioffs_orbslao)
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
      end
