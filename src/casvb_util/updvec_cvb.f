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
      subroutine updvec_cvb(upd,iorb,jorb,niprev,iprev,
     >  orbs,north,corth)
c  Find update for IORB as projection of JORB on allowed space
      implicit real*8 (a-h,o-z)
#include "ext_cvb.fh"
#include "main_cvb.fh"
#include "optze_cvb.fh"
#include "files_cvb.fh"
#include "print_cvb.fh"

#include "malloc_cvb.fh"
      dimension upd(norb)
      dimension iprev(niprev),orbs(norb,norb)
      dimension north(norb),corth(norb,niorth)
      dimension dum(1)

      i1 = mstackr_cvb(norb*norb)
      noffort=0
      do 100 ior=1,iorb-1
100   noffort=noffort+north(ior)
c  Collect all constraints and find span :
      call span0_cvb(norb,norb)
      if(north(iorb).gt.0)
     >  call span1_cvb(corth(1,1+noffort),north(iorb),dum,norb,0)
      do 200 i=1,niprev
200   call span1_cvb(orbs(1,iprev(i)),1,dum,norb,0)
      call span1_cvb(orbs(1,iorb),1,dum,norb,0)
      call span2_cvb(w(i1),ncon,dum,norb,0)

c  Orthogonalise update to all remaining constraints
      call fmove_cvb(orbs(1,jorb),upd,norb)
      call schmidtd_cvb(w(i1),ncon,upd,1,dum,norb,0)
      call mfreer_cvb(i1)
      return
      end
