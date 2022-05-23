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
      subroutine formats_cvb()
      implicit real*8 (a-h,o-z)
#include "main_cvb.fh"
#include "optze_cvb.fh"
#include "files_cvb.fh"
#include "print_cvb.fh"

#include "formats_cvb.fh"

      lenfld=4+iprec
      ix2=max(0,(lenfld-3)/2)
      ix1=max(0,lenfld-3-ix2)
c  iprec=8 => formMXP1='(4x,8(5x,i3,4x))'
      formMXP1='(4x,8('
      if(ix1.gt.0)then
        call appendint_cvb(formMXP1,ix1,0)
        call appendchr_cvb(formMXP1,'x,',0)
      endif
      call appendchr_cvb(formMXP1,'i3',0)
      if(ix2.gt.0)then
        call appendchr_cvb(formMXP1,',',0)
        call appendint_cvb(formMXP1,ix2,0)
        call appendchr_cvb(formMXP1,'x',0)
      endif
      call appendchr_cvb(formMXP1,'))',0)
c  iprec=8 => formMXP2='(4x,8(7x,i3,6x))'
      formMXP2='(4x,8('
      if(ix1.gt.0)then
        call appendint_cvb(formMXP2,ix1+2,0)
        call appendchr_cvb(formMXP2,'x,',0)
      endif
      call appendchr_cvb(formMXP2,'i3',0)
      if(ix2.gt.0)then
        call appendchr_cvb(formMXP2,',',0)
        call appendint_cvb(formMXP2,ix2+2,0)
        call appendchr_cvb(formMXP2,'x',0)
      endif
      call appendchr_cvb(formMXP2,'))',0)
c  iprec=8 => formMXP3='(1x,i3,8f12.8)'
      formMXP3='(1x,i3,8f'
      call appendint_cvb(formMXP3,lenfld,0)
      call appendchr_cvb(formMXP3,'.',0)
      call appendint_cvb(formMXP3,iprec,0)
      call appendchr_cvb(formMXP3,')',0)
c  iprec=8 => formMXP4='(1x,i3,8d16.8)'
      formMXP4='(1x,i3,8d'
      call appendint_cvb(formMXP4,lenfld+4,0)
      call appendchr_cvb(formMXP4,'.',0)
      call appendint_cvb(formMXP4,iprec,0)
      call appendchr_cvb(formMXP4,')',0)
c  iprec=8 => formMXP5='(4x,8f12.8)'
      formMXP5='(4x,8f'
      call appendint_cvb(formMXP5,lenfld,0)
      call appendchr_cvb(formMXP5,'.',0)
      call appendint_cvb(formMXP5,iprec,0)
      call appendchr_cvb(formMXP5,')',0)
c  iprec=8 => formMXP6='(4x,8d16.8)'
      formMXP6='(4x,8d'
      call appendint_cvb(formMXP6,lenfld+4,0)
      call appendchr_cvb(formMXP6,'.',0)
      call appendint_cvb(formMXP6,iprec,0)
      call appendchr_cvb(formMXP6,')',0)

c  iprec=8 => formE='(a,f16.10)'
      formE='(a,f'
      call appendint_cvb(formE,iprec+8,0)
      call appendchr_cvb(formE,'.',0)
      call appendint_cvb(formE,iprec+2,0)
      call appendchr_cvb(formE,')',0)
c  iprec=8 => formSymW='(a,i2,a,4e16.8)'
      formSymW='(a,i2,a,4e'
      call appendint_cvb(formSymW,iprec+8,0)
      call appendchr_cvb(formSymW,'.',0)
      call appendint_cvb(formSymW,iprec,0)
      call appendchr_cvb(formSymW,')',0)
c  iprec=8 => formVBWnorm='(a,f12.8,a)'
      formVBWnorm='(a,f'
      call appendint_cvb(formVBWnorm,iprec+4,0)
      call appendchr_cvb(formVBWnorm,'.',0)
      call appendint_cvb(formVBWnorm,iprec,0)
      call appendchr_cvb(formVBWnorm,',a)',0)

c  iprec=8 => formChk1='(a,d16.8)'
      formChk1='(a,d'
      call appendint_cvb(formChk1,iprec+8,0)
      call appendchr_cvb(formChk1,'.',0)
      call appendint_cvb(formChk1,iprec,0)
      call appendchr_cvb(formChk1,')',0)
c  iprec=8 => formChk2='(5(6x,a,1x))'
      lenfld=8+iprec
      ix1=max(0,min(lenfld-9,6))
      ix2=max(0,lenfld-9-ix1)
      formChk2='(5('
      call appendint_cvb(formChk2,ix1,0)
      call appendchr_cvb(formChk2,'x,a,',0)
      call appendint_cvb(formChk2,ix2,0)
      call appendchr_cvb(formChk2,'x))',0)
c  iprec=8 => formChk3='(5d16.8)'
      formChk3='(5d'
      call appendint_cvb(formChk3,iprec+8,0)
      call appendchr_cvb(formChk3,'.',0)
      call appendint_cvb(formChk3,iprec,0)
      call appendchr_cvb(formChk3,')',0)
c  iprec=8 => formcvp='(2(a,d16.8))'
      formcvp='(2(a,d'
      call appendint_cvb(formcvp,iprec+8,0)
      call appendchr_cvb(formcvp,'.',0)
      call appendint_cvb(formcvp,iprec,0)
      call appendchr_cvb(formcvp,'))',0)
c  iprec=8 => formAD='(a,10d16.8)'
      formAD='(a,10d'
      call appendint_cvb(formAD,iprec+8,0)
      call appendchr_cvb(formAD,'.',0)
      call appendint_cvb(formAD,iprec,0)
      call appendchr_cvb(formAD,')',0)
c  iprec=8 => formAF='(a,10f16.8)'
      formAF='(a,10f'
      call appendint_cvb(formAF,iprec+8,0)
      call appendchr_cvb(formAF,'.',0)
      call appendint_cvb(formAF,iprec,0)
      call appendchr_cvb(formAF,')',0)
c  iprec=8 => form2AD='(2(a,d16.8))'
      form2AD='(2(a,d'
      call appendint_cvb(form2AD,iprec+8,0)
      call appendchr_cvb(form2AD,'.',0)
      call appendint_cvb(form2AD,iprec,0)
      call appendchr_cvb(form2AD,'))',0)
c  iprec=8 => form2AF='(2(a,f16.8))'
      form2AF='(2(a,f'
      call appendint_cvb(form2AF,iprec+8,0)
      call appendchr_cvb(form2AF,'.',0)
      call appendint_cvb(form2AF,iprec,0)
      call appendchr_cvb(form2AF,'))',0)
c  iprec=8 => formroot='(/,a,i3,a,f16.10,a)'
      formroot='(/,a,i3,a,f'
      call appendint_cvb(formroot,iprec+8,0)
      call appendchr_cvb(formroot,'.',0)
      call appendint_cvb(formroot,iprec+2,0)
      call appendchr_cvb(formroot,',a)',0)
      return
      end
