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
c  Contents of CI vectors :
      subroutine setcnt_cvb(xident_ci,idep1)
      implicit real*8 (a-h,o-z)
#include "ext_cvb.fh"
#include "main_cvb.fh"
#include "optze_cvb.fh"
#include "files_cvb.fh"
#include "print_cvb.fh"
      dimension xident_ci(1)

      ident_ci=nint(xident_ci(1))
      call setcnt2_cvb(ident_ci,idep1)
      return
      end
      subroutine setcnt2_cvb(ident_ci,idep1)
      implicit real*8 (a-h,o-z)
#include "ext_cvb.fh"
#include "main_cvb.fh"
#include "optze_cvb.fh"
#include "files_cvb.fh"
#include "print_cvb.fh"

      icnt_ci(ident_ci)=idep1
      return
      end
      function igetcnt2_cvb(ident_ci)
      implicit real*8 (a-h,o-z)
#include "ext_cvb.fh"
#include "main_cvb.fh"
#include "optze_cvb.fh"
#include "files_cvb.fh"
#include "print_cvb.fh"

      igetcnt2_cvb=icnt_ci(ident_ci)
      return
      end
      logical function tstcnt_cvb(xident_ci,idep1)
      implicit real*8 (a-h,o-z)
#include "main_cvb.fh"
#include "optze_cvb.fh"
#include "files_cvb.fh"
#include "print_cvb.fh"
      dimension xident_ci(1)

      ident_ci=nint(xident_ci(1))
      tstcnt_cvb=icnt_ci(ident_ci).eq.idep1
      return
      end
