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
      subroutine casinfo2_cvb()
      implicit real*8 (a-h,o-z)
#include "ext_cvb.fh"
#include "main_cvb.fh"
#include "optze_cvb.fh"
#include "files_cvb.fh"
#include "print_cvb.fh"

#include "casinfo_cvb.fh"

c  Information from molcas interface file "JOBIPH" :
      call rdjobiph_cvb('JOBIPH')
      call setjobiph_cvb(iorcore_c,iorclos_c,iorocc_c,mxirrep,
     >  nstsym_c,weight_c,istnel_c,istsy_c,istms2_c,nstats_c,
     >  mxstt_ci,mxstsy_ci,nel_c,norb_c,i2s_c,isym_c,mcore_c,neltot_c)
      return
      end
