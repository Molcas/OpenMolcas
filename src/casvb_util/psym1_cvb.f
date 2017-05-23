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
      subroutine psym1_cvb(civec1,civec2,osym,ientry)
      implicit real*8 (a-h,o-z)
#include "ext_cvb.fh"
#include "main_cvb.fh"
#include "optze_cvb.fh"
#include "files_cvb.fh"
#include "print_cvb.fh"

#include "malloc_cvb.fh"
      dimension civec1(nda,ndb),civec2(nda,ndb),osym(mxirrep)

      k1 = mstacki_cvb(nda)
      k2 = mstacki_cvb(ndb)
      k3 = mstacki_cvb(mxirrep+1)
      k4 = mstacki_cvb(mxirrep+1)
      k5 = mstacki_cvb(nda)
      k6 = mstacki_cvb(ndb)
      k7 = mstacki_cvb(mxirrep)
      k8 = mstacki_cvb(mxirrep)
      k9 = mstacki_cvb(mxirrep)
      k10= mstacki_cvb(norb+1)
      k11= mstacki_cvb(norb+1)
      k12= mstacki_cvb(norb+1)
      k13= mstacki_cvb(norb+1)
      k14= mstacki_cvb(norb+1)
      k15= mstacki_cvb((norb+1)*(nalf+1))
      k16= mstacki_cvb((norb+1)*(nbet+1))
      k17= mstacki_cvb(mxirrep)

      call symgen_cvb(nalf,nbet,nda,ndb,
     >  iw(k1),iw(k2),iw(k3),iw(k4),
     >  iw(k5),iw(k6),iw(k7),iw(k8),iw(k9),
     >  iw(k10),iw(k11),iw(k12),iw(k13),iw(k14),iw(k15),
     >  iw(k16),iw(k17))

      call mfreei_cvb(k5)
      call psym2_cvb(civec1,civec2,
     >  iw(k1),iw(k2),iw(k3),iw(k4),osym,ientry)
      call mfreei_cvb(k1)
      return
      end
