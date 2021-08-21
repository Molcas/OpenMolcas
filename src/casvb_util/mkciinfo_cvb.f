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
      subroutine mkciinfo_cvb()
      implicit real*8 (a-h,o-z)
#include "main_cvb.fh"
#include "optze_cvb.fh"
#include "files_cvb.fh"
#include "print_cvb.fh"

#include "WrkSpc.fh"

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
     >  iwork(ll(1)),iwork(ll(2)),iwork(ll(3)),iwork(ll(4)),
     >  iwork(ll(5)),iwork(ll(6)),
     >  work(ll(9)),work(ll(10)),
     >  iwork(k1),iwork(k2),iwork(k3),iwork(k4),iwork(k5),iwork(k6),
     >  iwork(k7),iwork(k8),iwork(k9),iwork(k10),iwork(k11))
      call mfreei_cvb(k1)
      return
      end
