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
      subroutine evb2cas_cvb(orbs,cvb,fx,ioptc,iter)
      implicit real*8 (a-h,o-z)
c ... Files/Hamiltonian available ...
      logical, external :: tstfile_cvb
#include "main_cvb.fh"
#include "optze_cvb.fh"
#include "files_cvb.fh"
#include "print_cvb.fh"

#include "WrkSpc.fh"
      dimension orbs(norb,norb),cvb(nvb)

      if(tstfile_cvb(66000.2d0))then
        i1=mstackr_cvb(norb*norb+nvb)
        call rdr_cvb(work(i1),norb*norb+nvb,66000.2d0,0)
        call subvec(work(i1),orbs,work(i1),norb*norb)
        call subvec(work(norb*norb+i1),cvb,work(norb*norb+i1),nvb)
        dxnrm=dnrm2_(norb*norb+nvb,work(i1),1)
        call findamx_cvb(work(i1),norb*norb+nvb,dx_amx,idum)
        call mfreer_cvb(i1)
      endif
      call wrr_cvb(orbs,norb*norb,66000.2d0,0)
      call wrr_cvb(cvb,nvb,66000.2d0,norb*norb)
      call evb2cas2_cvb(orbs,cvb,ioptc,iter,fx,
     >  dxnrm,dx_amx,
     >  work(lc(1)),work(lc(2)),work(lc(3)),work(lc(4)),work(lc(5)),
     >  work(lw(2)),work(lw(3)))
      return
      end
