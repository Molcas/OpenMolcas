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
c  *  GETNCI    := Get number of CASSCF determinants in each irrep.    *
c  *                                                                   *
c  *********************************************************************
      subroutine getnci_cvb(nciloc,nelloc,i2sloc,isymloc)
      implicit real*8 (a-h,o-z)
#include "main_cvb.fh"
#include "optze_cvb.fh"
#include "files_cvb.fh"
#include "print_cvb.fh"

#include "WrkSpc.fh"
      dimension nciloc(*)

      nalf1=(nelloc+i2sloc)/2
      nbet1=nelloc-nalf1
      call icomb_cvb(norb,nalf1,nda1)
      call icomb_cvb(norb,nbet1,ndb1)
      k1 = mstacki_cvb(nda1)
      k2 = mstacki_cvb(ndb1)
      k3 = mstacki_cvb(mxirrep+1)
      k4 = mstacki_cvb(mxirrep+1)
      k5 = mstacki_cvb(nda1)
      k6 = mstacki_cvb(ndb1)
      k7 = mstacki_cvb(mxirrep)
      k8 = mstacki_cvb(mxirrep)
      k9 = mstacki_cvb(mxirrep)
      k10= mstacki_cvb(norb+1)
      k11= mstacki_cvb(norb+1)
      k12= mstacki_cvb(norb+1)
      k13= mstacki_cvb(norb+1)
      k14= mstacki_cvb(norb+1)
      k15= mstacki_cvb((norb+1)*(nalf1+1))
      k16= mstacki_cvb((norb+1)*(nbet1+1))
      k17= mstacki_cvb(mxirrep)

      call symgen_cvb(nalf1,nbet1,nda1,ndb1,
     >  iwork(k1),iwork(k2),iwork(k3),iwork(k4),
     >  iwork(k5),iwork(k6),iwork(k7),iwork(k8),iwork(k9),
     >  iwork(k10),iwork(k11),iwork(k12),iwork(k13),iwork(k14),
     >  iwork(k15),
     >  iwork(k16),iwork(k17))

      if(isymloc.eq.0)then
        call imove_cvb(iwork(k7),nciloc,mxirrep)
      else
        nciloc(1)=iwork(isymloc+k7-1)
      endif
      call mfreei_cvb(k1)
      return
      end
