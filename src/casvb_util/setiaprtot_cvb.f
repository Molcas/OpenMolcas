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
      subroutine setiaprtot_cvb()
      implicit real*8 (a-h,o-z)
#include "ext_cvb.fh"
#include "main_cvb.fh"
#include "optze_cvb.fh"
#include "files_cvb.fh"
#include "print_cvb.fh"


#include "frag_cvb.fh"
#include "malloc_cvb.fh"
      dimension dum1(1), dum2(1)

      k1=mstackr_cvb(nda*ndb)
      call dpci2vb_cvb(w(k1),dum1,dum2,0,dum3,4)
      call setiaprtot2_cvb(w(k1),
     >  iw(ll(11)),iw(ll(12)),iw(ll(13)),iw(ll(14)),npvb,
     >  nda,ndb)
      call mfreer_cvb(k1)
      return
      end
