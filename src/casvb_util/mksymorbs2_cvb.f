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
      subroutine mksymorbs2_cvb(orbs,sorbs)
      implicit real*8 (a-h,o-z)
#include "main_cvb.fh"
#include "optze_cvb.fh"
#include "files_cvb.fh"
#include "print_cvb.fh"

      dimension orbs(norb,norb)
      dimension sorbs(norb,norb)
      dimension dum(1)
      save thresh
      data thresh/1.d-7/

      if(sym)then
        call fmove_cvb(orbs,sorbs,norb*norb)
        nconstr_kp=nconstr
        nconstr=0
        call symtrizorbs_cvb(orbs)
        nconstr=nconstr_kp
        call subvec(sorbs,orbs,sorbs,norb*norb)
        delorbs=dnrm2_(norb*norb,sorbs,1)
        if(delorbs.gt.thresh.and.ip(1).ge.2)then
          write(6,'(/,a)') ' Change in symmetrized orbitals:'
          call report_cvb(sorbs,norb)
        endif
        call nize_cvb(orbs,norb,dum,norb,0,0)
        if(delorbs.gt.thresh.and.ip(1).ge.2)then
          write(6,'(a)') ' Orbitals after symmetrization:'
          call report_cvb(orbs,norb)
        endif
        if(abs(detm_cvb(orbs,norb)).lt.1d-8)then
          write(6,*)' Fatal error - orbital matrix singular',
     >      ' after symmetrization!'
          call abend_cvb()
        endif
      endif
      return
      end
