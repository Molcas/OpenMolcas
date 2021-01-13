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
c  ************************************
c  ** Print and save final VB energy **
c  ************************************
      subroutine finalresult_cvb()
      implicit real*8 (a-h,o-z)
c ... Make: up to date? ...
      logical, external :: up2date_cvb
#include "main_cvb.fh"
#include "optze_cvb.fh"
#include "files_cvb.fh"
#include "print_cvb.fh"

      common /nfin_comcvb/nfinal

      nfinal=nfinal+1

      if((.not.variat).and.up2date_cvb('SVB'))then
        call add_info('SVB',[abs(svb)],1,7)
      endif
      if((.not.variat).and.up2date_cvb('EVB'))then
        call add_info('EVB',[evb],1,7)
      endif
      return
      entry finalresult_init_cvb()
      nfinal=0
      return
      end
