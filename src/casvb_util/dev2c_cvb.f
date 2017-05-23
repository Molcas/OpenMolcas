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
      subroutine dev2c_cvb(v1,cfrom,hessorb,oaa2)
c  Calculate V1 EijEkl CFROM
      implicit real*8 (a-h,o-z)
#include "ext_cvb.fh"
#include "main_cvb.fh"
#include "optze_cvb.fh"
#include "files_cvb.fh"
#include "print_cvb.fh"

#include "malloc_cvb.fh"
      dimension hessorb(nprorb,nprorb),v1(*),cfrom(*)

      iv1=nint(v1(1))
      icfrom=nint(cfrom(1))
      n_2el=n_2el+1
      if(iform_ci(icfrom).ne.0)then
        write(6,*)' Unsupported format in DEV2C :',iform_ci(icfrom)
        call abend_cvb()
      endif

      call dev2c_2_cvb(w(iaddr_ci(iv1)),w(iaddr_ci(icfrom)),
     >  hessorb,oaa2,nprorb,
     >  iw(ll(1)),iw(ll(2)),iw(ll(3)),iw(ll(4)),iw(ll(5)),iw(ll(6)),
     >  w(ll(9)),w(ll(10)),
     >  iw(ll(11)),iw(ll(12)),iw(ll(13)),iw(ll(14)),npvb,
     >  nda,ndb,n1a,n1b,nam1,nbm1,norb,projcas,sc,absym(3))
      return
      end
