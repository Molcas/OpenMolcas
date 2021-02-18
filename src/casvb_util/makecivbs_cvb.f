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
      subroutine makecivbs_cvb(civbs,orbs,gjorb,gjorb2,gjorb3,cvbdet)
c  Construct CIVBS ( = T(s) * CIVB ) :
      implicit real*8 (a-h,o-z)
c ... Content of CI vectors ...
      logical, external :: tstcnt_cvb
#include "main_cvb.fh"
#include "optze_cvb.fh"
#include "files_cvb.fh"
#include "print_cvb.fh"

      dimension orbs(norb,norb)
      dimension civbs(ndet)
      dimension gjorb(*),gjorb2(*),gjorb3(*)
      dimension cvbdet(ndetvb)

      if(tstcnt_cvb(civbs,4))return

      call vb2cic_cvb(cvbdet,civbs)
      call applyts_cvb(civbs,orbs,gjorb,gjorb2,gjorb3)
      call setcnt_cvb(civbs,4)
      return
      end
