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
c  *********************************************************************
c  *                                                                   *
c  *  DEV2A  := calculate two-electron Hessian                         *
c  *                                                                   *
c  *********************************************************************
      subroutine dev2a_cvb(v1,v2,cfrom,hessorb,oaa2,aa1)
c  Calculates V1 EijEkl CFROM and V2 EijEkl CFROM
      implicit real*8 (a-h,o-z)
#include "main_cvb.fh"
#include "optze_cvb.fh"
#include "files_cvb.fh"
#include "print_cvb.fh"

#include "WrkSpc.fh"
      dimension hessorb(nprorb,nprorb),v1(*),v2(*),cfrom(*)

      iv1=nint(v1(1))
      iv2=nint(v2(1))
      icfrom=nint(cfrom(1))
      n_2el=n_2el+2
      if(iform_ci(icfrom).ne.0)then
        write(6,*)' Unsupported format in DEV2A :',iform_ci(icfrom)
        call abend_cvb()
      endif

      call dev2a_2_cvb(work(iaddr_ci(iv1)),work(iaddr_ci(iv2)),
     >  work(iaddr_ci(icfrom)),
     >  hessorb,oaa2,aa1,nprorb,
     >  iwork(ll(1)),iwork(ll(2)),iwork(ll(3)),iwork(ll(4)),
     >  iwork(ll(5)),iwork(ll(6)),
     >  work(ll(9)),work(ll(10)),
     >  iwork(ll(11)),iwork(ll(12)),iwork(ll(13)),iwork(ll(14)),npvb,
     >  nda,ndb,n1a,n1b,nam1,nbm1,norb,projcas,sc,absym(3))
      return
      end
