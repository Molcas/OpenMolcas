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
      subroutine casinfodef_cvb()
      implicit real*8 (a-h,o-z)
#include "ext_cvb.fh"
#include "main_cvb.fh"
#include "optze_cvb.fh"
#include "files_cvb.fh"
#include "print_cvb.fh"

#include "casinfo_cvb.fh"
#include "inpmod_cvb.fh"

c  Counters
      nstsym_d=0

      noe=2*mxorb
      if(.not.variat)then
        strtci=recn_jobold
      else
        strtci=recn_jobiph
      endif
      strtmo=recn_jobiph
      strtint=recn_oneint
      strtvb=recn_vbwfn
      savvb=recn_vbwfn
      savvbci=recn_jobiph

      if(inputmode.eq.2)then
        do 100 i=1,mxstsy_ci
        iorcore_d(i)=-1
        iorclos_d(i)=-1
100     iorocc_d(i)=-1
      endif
      return
      end
