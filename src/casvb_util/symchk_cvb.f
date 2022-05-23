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
      subroutine symchk_cvb()
      implicit real*8 (a-h,o-z)
c ... Make: up to date? ...
      logical, external :: up2date_cvb
#include "main_cvb.fh"
#include "optze_cvb.fh"
#include "files_cvb.fh"
#include "print_cvb.fh"

      logical recinpcmp_cvb

      if(up2date_cvb('SYMINIT'))then
c  iorts:
        if(recinpcmp_cvb(14))call touch_cvb('ORBFREE')
c  irots:
        if(recinpcmp_cvb(15))call touch_cvb('ORBFREE')
c  iorbrel:
        if(recinpcmp_cvb(9))then
          call touch_cvb('SYMINIT')
          call touch_cvb('ORBFREE')
        endif
c  ifxorb:
        if(recinpcmp_cvb(11))then
          call touch_cvb('SYMINIT')
          call touch_cvb('ORBFREE')
        endif
      endif
      if(up2date_cvb('CONSTRUC'))then
c  ifxstr:
        if(recinpcmp_cvb(12))then
          call touch_cvb('CONSTRUC')
          call touch_cvb('CIFREE')
        endif
c  idelstr:
        if(recinpcmp_cvb(13))then
          call touch_cvb('CONSTRUC')
          call touch_cvb('CIFREE')
        endif
c  izeta:
        if(recinpcmp_cvb(16))then
          call touch_cvb('CONSTRUC')
          call touch_cvb('CIFREE')
        endif
      endif
      return
      end
