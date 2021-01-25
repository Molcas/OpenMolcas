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
c  ****************************************************
c  ** Routines that deal with spatial configurations **
c  ****************************************************
c  *********************************************************************
c  *                                                                   *
c  *  CNFCHECK := Check input configurations and convert to occ no     *
c  *              representation if necessary.                         *
c  *              Also sort in order of increasing ionicity.           *
c  *                                                                   *
c  *********************************************************************
      subroutine cnfcheck_cvb(iconfs,nconf1,nel1)
      implicit real*8 (a-h,o-z)
#include "main_cvb.fh"

#include "optze_cvb.fh"
#include "files_cvb.fh"
#include "print_cvb.fh"


#include "malloc_cvb.fh"
      dimension iconfs(noe,nconf1)

      i1 = mstacki_cvb(noe)
      call cnfcheck2_cvb(iconfs,nconf1,nel1,iw(i1))
      call mfreei_cvb(i1)
      ioncty = mstacki_cvb(nconf1)
      iconfs2 = mstacki_cvb(noe*nconf1)
      call cnfsort_cvb(iconfs,nconf1,nel1,iw(ioncty),iw(iconfs2))
      call mfreei_cvb(ioncty)
      return
      end
