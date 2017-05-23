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
      subroutine wrgspr_cvb(recn,c,i,n,ic,ierr)
      implicit real*8 (a-h,o-z)
#include "ext_cvb.fh"
#include "main_cvb.fh"
#include "optze_cvb.fh"
#include "files_cvb.fh"
#include "print_cvb.fh"

#include "malloc_cvb.fh"
#include "mo_cvb.fh"
      dimension c(n)

c  Read header :
      call rdheader_cvb(recn,norb1,nbas_mo1,nvb1,kbasiscvb1,
     >  ioffs_orbs,ioffs_cvb,ioffs_orbsao,ioffs_orbslao)

      if(ic.eq.1)then
c  >>> Orbital write, I is orbital number :
        if(i.gt.norb1)then
          ierr=1
          return
        endif
        ioffs=(i-1)*norb1+ioffs_orbs
        call wrr_cvb(c,min(norb1,n),recn,ioffs)
      elseif(ic.eq.2)then
c  >>> Structure write, I is starting structure coefficient :
        if(i.gt.nvb1)then
          ierr=1
          return
        endif
        ioffs=i-1+ioffs_cvb
        call wrr_cvb(c,min(nvb1,n),recn,ioffs)
      elseif(ic.eq.3)then
c  >>> Write of orbital in AO basis, I is orbital number :
        if(i.gt.norb1)then
          ierr=1
          return
        endif
        if(nbas_mo1.eq.0)then
          nbas_mo1=nbas_mo
          call wrheader_cvb(recn,norb1,nbas_mo1,nvb1,kbasiscvb1,
     >      ioffs_orbs,ioffs_cvb,ioffs_orbsao,ioffs_orbslao)
        endif
c  Error return if AO bases are not identical :
        if(nbas_mo1.ne.n)then
          ierr=1
          return
        endif
        ioffs=(i-1)*nbas_mo1+ioffs_orbsao
        call wrr_cvb(c,min(nbas_mo1,n),recn,ioffs)
      elseif(ic.eq.4)then
c  >>> Write of localized orbital in AO basis, I is orbital number :
        if(i.gt.norb1)then
          ierr=1
          return
        endif
        if(nbas_mo1.eq.0)then
          nbas_mo1=nbas_mo
          call wrheader_cvb(recn,norb1,nbas_mo1,nvb1,kbasiscvb1,
     >      ioffs_orbs,ioffs_cvb,ioffs_orbsao,ioffs_orbslao)
        endif
c  Error return if AO bases are not identical :
        if(nbas_mo1.ne.n)then
          ierr=1
          return
        endif
        ioffs=(i-1)*nbas_mo1+ioffs_orbslao
        call wrr_cvb(c,min(nbas_mo1,n),recn,ioffs)
      endif
      return
      end
