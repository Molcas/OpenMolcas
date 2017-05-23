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
      subroutine mkciinfo_cvb()
      implicit real*8 (a-h,o-z)
#include "ext_cvb.fh"
#include "main_cvb.fh"
#include "optze_cvb.fh"
#include "files_cvb.fh"
#include "print_cvb.fh"

#include "malloc_cvb.fh"

      k1 = mstacki_cvb((norb+1)*(nalf+1))
      k2 = mstacki_cvb((norb+1)*(nbet+1))
      k3 = mstacki_cvb((norb+1)*nalf)
      k4 = mstacki_cvb((norb+1)*nbet)
      k5 = mstacki_cvb(norb+1)
      k6 = mstacki_cvb(norb+1)
      k7 = mstacki_cvb(norb+1)
      k8 = mstacki_cvb(norb+1)
      k9 = mstacki_cvb(norb+1)
      k10= mstacki_cvb(norb)
      k11= mstacki_cvb(norb)
      call mkciinfo2_cvb(
     >  iw(ll(1)),iw(ll(2)),iw(ll(3)),iw(ll(4)),iw(ll(5)),iw(ll(6)),
     >  w(ll(9)),w(ll(10)),
     >  iw(k1),iw(k2),iw(k3),iw(k4),iw(k5),iw(k6),
     >  iw(k7),iw(k8),iw(k9),iw(k10),iw(k11))
      call mfreei_cvb(k1)
      return
      end
