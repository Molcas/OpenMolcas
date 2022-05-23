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
      subroutine symtrizcvb2_cvb(vecstr,
     >  izeta,ipermzeta,dvbdet,vecstr2)
      implicit real*8 (a-h,o-z)
#include "main_cvb.fh"
#include "optze_cvb.fh"
#include "files_cvb.fh"
#include "print_cvb.fh"

      dimension vecstr(nvb)
      dimension izeta(nsyme),ipermzeta(norb,nzeta)
      dimension dvbdet(ndetvb),vecstr2(nvb)

      izeta1=0
      do 100 isyme=1,nsyme
      if(izeta(isyme).ne.0)then
        izeta1=izeta1+1
        call str2vbc_cvb(vecstr,dvbdet)
        call permvb_cvb(dvbdet,ipermzeta(1,izeta1))
        call vb2strc_cvb(dvbdet,vecstr2)
        call daxpy_(nvb,dble(izeta(isyme)),vecstr2,1,vecstr,1)
      endif
100   continue
      if(izeta1.gt.0)call dscal_(nvb,1d0/DBLE(2**izeta1),vecstr,1)
      return
      end
