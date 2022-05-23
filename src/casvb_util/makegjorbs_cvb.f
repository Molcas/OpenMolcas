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
      subroutine makegjorbs_cvb(orbs,gjorb,gjorb2,gjorb3)
c  Construct Gauss-Jordan factorizations of ORBS, ORBS transpose,
c  and overlap matrix corresonding to ORBS :
      implicit real*8 (a-h,o-z)
#include "main_cvb.fh"
#include "optze_cvb.fh"
#include "files_cvb.fh"
#include "print_cvb.fh"

#include "malloc_cvb.fh"
      dimension orbs(norb,norb),gjorb(*),gjorb2(*),gjorb3(*)

      iowrk  = mstackr_cvb(norb*norb)

      call gaussj_cvb(orbs,gjorb)

      call transp_cvb(orbs,w(iowrk),norb,norb)
      call gaussj_cvb(w(iowrk),gjorb2)

      call mxattb_cvb(orbs,orbs,norb,norb,norb,w(iowrk))
      call gaussj_cvb(w(iowrk),gjorb3)

      call mfreer_cvb(iowrk)

      return
      end
