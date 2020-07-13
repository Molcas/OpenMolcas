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
      subroutine orthcvb_cvb(c,nparm1)
      implicit real*8 (a-h,o-z)
#include "ext_cvb.fh"
#include "main_cvb.fh"
#include "optze_cvb.fh"
#include "files_cvb.fh"
#include "print_cvb.fh"

#include "frag_cvb.fh"
#include "malloc_cvb.fh"
      dimension c(*)

      ioffs=nparm1-nprvb+1
      if(nfrag.le.1)then
        call daxpy_(nprvb,-ddot_(nprvb,w(lv(2)),1,c(ioffs),1)/cvbnrm,
     >    w(lv(2)),1,c(ioffs),1)
      else
        ifr_off=0
        do 100 ifrag=1,nfrag
        call daxpy_(nvb_fr(ifrag),
     >    -ddot_(nvb_fr(ifrag),w(ifr_off+lv(2)),1,c(ifr_off+ioffs),1)
     >    /cvbnrm_fr(ifrag),
     >   w(ifr_off+lv(2)),1,c(ifr_off+ioffs),1)
        ifr_off=ifr_off+nvb_fr(ifrag)
100     continue
      endif
      return
      entry orthcvb_init_cvb()
      if(nfrag.le.1)then
        cvbnrm=ddot_(nvb,w(lv(2)),1,w(lv(2)),1)
      else
        ifr_off=0
        do 200 ifrag=1,nfrag
        cvbnrm_fr(ifrag)=ddot_(nvb_fr(ifrag),w(ifr_off+lv(2)),1,
     >    w(ifr_off+lv(2)),1)
        ifr_off=ifr_off+nvb_fr(ifrag)
200     continue
      endif
      return
      end
