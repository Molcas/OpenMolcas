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
      subroutine setstrtvb_cvb(recn)
      implicit real*8 (a-h,o-z)
c ... Files/Hamiltonian available ...
      logical, external :: tstfile_cvb
#include "main_cvb.fh"
#include "optze_cvb.fh"
#include "files_cvb.fh"
#include "print_cvb.fh"

#include "malloc_cvb.fh"
      save recdef
      data recdef/3200.2d0/

      if(recn.ne.zero)return
      if(.not.tstfile_cvb(recdef))return
      do 100 iadd=1,99
      if(.not.tstfile_cvb(recdef+DBLE(iadd)))then
        jadd=iadd-1
        recn=recdef+DBLE(jadd)
        return
      endif
100   continue
      return
      end
      subroutine setsavvb_cvb(recn)
      implicit real*8 (a-h,o-z)
c ... Files/Hamiltonian available ...
      logical, external :: tstfile_cvb
#include "main_cvb.fh"
#include "optze_cvb.fh"
#include "files_cvb.fh"
#include "print_cvb.fh"

#include "malloc_cvb.fh"
      save recdef
      data recdef/3200.2d0/
      if(recn.ne.zero)return
      do 200 iadd=0,99
      if(.not.tstfile_cvb(recdef+DBLE(iadd)))then
        recn=recdef+DBLE(iadd)
        return
      endif
200   continue
      recn=recdef+DBLE(99)
      return
      end
      logical function ifmos_cvb()
      implicit real*8 (a-h,o-z)

      ifmos_cvb=.true.
      return
      end
      logical function ifcasci_cvb()
      implicit real*8 (a-h,o-z)
#include "main_cvb.fh"
#include "optze_cvb.fh"
#include "files_cvb.fh"
#include "print_cvb.fh"

#include "applyh_cvb.fh"

      call f_inquire('JOBOLD',ifcasci_cvb)
c  In variational calculations, CI vectors will be of no use
c  unless they correspond to MOs identical to the present ones :
      if(variat.and.(invec_cvb.ne.3.or.nmcscf.gt.1))ifcasci_cvb=.false.
      return
      end
      logical function ifhamil_cvb()
      implicit real*8 (a-h,o-z)
      ifhamil_cvb=.true.
      return
      end
      logical function valid_cvb(fileid)
      implicit real*8 (a-h,o-z)

      valid_cvb=(fileid.ge.1d-2)
      return
      end
