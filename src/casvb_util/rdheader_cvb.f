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
      subroutine rdheader_cvb(recn,norb1,nbas_mo1,nvb1,kbasiscvb1,
     >  ioffs_orbs,ioffs_cvb,ioffs_orbsao,ioffs_orbslao)
      implicit real*8 (a-h,o-z)
#include "main_cvb.fh"
#include "optze_cvb.fh"
#include "files_cvb.fh"
#include "print_cvb.fh"

#include "malloc_cvb.fh"
      logical debug
      data debug/.false./

      iheader=mstackiz_cvb(10)
      call mfreei_cvb(iheader)
      iheader=mstackiz_cvb(10)
      call rdi_cvb(iw(iheader),10,recn,0)
      norb1=iw(iheader)
      nbas_mo1=iw(1+iheader)
      nvb1=iw(2+iheader)
      kbasiscvb1=iw(3+iheader)
      ioffs_orbs=iw(5+iheader)
      ioffs_cvb=iw(6+iheader)
      ioffs_orbsao=iw(7+iheader)
      ioffs_orbslao=iw(8+iheader)
      call mfreei_cvb(iheader)
      if(debug)then
        write(6,*)' rdheader :'
        write(6,*)' ----------'
        write(6,*)' norb1         :',norb1
        write(6,*)' nbas_mo1      :',nbas_mo1
        write(6,*)' nvb1          :',nvb1
        write(6,*)' kbasiscvb1    :',kbasiscvb1
        write(6,*)' ioffs_orbs    :',ioffs_orbs
        write(6,*)' ioffs_cvb     :',ioffs_cvb
        write(6,*)' ioffs_orbsao  :',ioffs_orbsao
        write(6,*)' ioffs_orbslao :',ioffs_orbslao
      endif
      return
      end
      subroutine wrheader_cvb(recn,norb1,nbas_mo1,nvb1,kbasiscvb1,
     >  ioffs_orbs,ioffs_cvb,ioffs_orbsao,ioffs_orbslao)
      implicit real*8 (a-h,o-z)
#include "main_cvb.fh"
#include "optze_cvb.fh"
#include "files_cvb.fh"
#include "print_cvb.fh"

#include "malloc_cvb.fh"
      logical debug
      data debug/.false./

        ioffs_orbs=ihlf_cvb(10)
        ioffs_cvb=ioffs_orbs+norb1*norb1
        ioffs_orbsao=ioffs_cvb+nvb1
        ioffs_orbslao=ioffs_orbsao+norb1*nbas_mo1
        call reserv_cvb(ioffs_orbslao+norb1*nbas_mo1,recn)
      iheader=mstackiz_cvb(10)
      iw(iheader)=norb1
      iw(1+iheader)=nbas_mo1
      iw(2+iheader)=nvb1
      iw(3+iheader)=kbasiscvb1
      iw(5+iheader)=ioffs_orbs
      iw(6+iheader)=ioffs_cvb
      iw(7+iheader)=ioffs_orbsao
      iw(8+iheader)=ioffs_orbslao
      call wri_cvb(iw(iheader),10,recn,0)
      call mfreei_cvb(iheader)
      if(debug)then
        write(6,*)' wrheader :'
        write(6,*)' ----------'
        write(6,*)' norb1         :',norb1
        write(6,*)' nbas_mo1      :',nbas_mo1
        write(6,*)' nvb1          :',nvb1
        write(6,*)' kbasiscvb1    :',kbasiscvb1
        write(6,*)' ioffs_orbs    :',ioffs_orbs
        write(6,*)' ioffs_cvb     :',ioffs_cvb
        write(6,*)' ioffs_orbsao  :',ioffs_orbsao
        write(6,*)' ioffs_orbslao :',ioffs_orbslao
      endif
      return
      end
