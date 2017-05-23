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
      subroutine symweight_cvb(civec1,civec2,osym)
      implicit real*8 (a-h,o-z)
#include "ext_cvb.fh"
#include "main_cvb.fh"
#include "optze_cvb.fh"
#include "files_cvb.fh"
#include "print_cvb.fh"

#include "malloc_cvb.fh"
      dimension osym(mxirrep),civec1(*),civec2(*)
c  *********************************************************************
c  *                                                                   *
c  *  SYMWEIGHT := CASSCF scalar product divided into irreps.          *
c  *                                                                   *
c  *********************************************************************

      icivec1=nint(civec1(1))
      icivec2=nint(civec2(1))
      call psym1_cvb(w(iaddr_ci(icivec1)),w(iaddr_ci(icivec2)),osym,2)
      return
      end
