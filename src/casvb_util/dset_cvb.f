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
      subroutine dset_cvb(iorbrel,ifxorb,ifxstr,
     >  idelstr,iorts,irots,izeta)
      implicit real*8 (a-h,o-z)
#include "ext_cvb.fh"
#include "main_cvb.fh"
#include "optze_cvb.fh"
#include "files_cvb.fh"
#include "print_cvb.fh"

#include "inpmod_cvb.fh"
#include "malloc_cvb.fh"
      dimension iorbrel(ndimrel),ifxorb(mxorb),
     >  iorts(*),irots(*),izeta(*)

c  Check if any molecular interaction constraints :
      if(ploc)call plcconst_plc()
      sym=(norbrel.gt.0.or.nort.gt.0.or.plc_const)
      call rdioff_cvb(9,recinp,ioffs)
      call wris_cvb(iorbrel,ndimrel,recinp,ioffs)
      call wrioff_cvb(10,recinp,ioffs)
      call wrioff_cvb(11,recinp,ioffs)
      call wris_cvb(ifxorb,norb,recinp,ioffs)
      call wrioff_cvb(12,recinp,ioffs)
      call wris_cvb(iw(ifxstr),nfxvb,recinp,ioffs)
      call wrioff_cvb(13,recinp,ioffs)
      call wris_cvb(iw(idelstr),nzrvb,recinp,ioffs)
      call wrioff_cvb(14,recinp,ioffs)
      call wris_cvb(iorts,2*nort,recinp,ioffs)
      call wrioff_cvb(15,recinp,ioffs)
      call wris_cvb(irots,2*ndrot,recinp,ioffs)
      call wrioff_cvb(16,recinp,ioffs)
      call wris_cvb(izeta,nsyme,recinp,ioffs)
      call wrioff_cvb(17,recinp,ioffs)
      return
      end
