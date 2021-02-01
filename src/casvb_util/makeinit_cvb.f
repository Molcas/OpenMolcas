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
c  *************************************
c  ** Routines to emulate unix "make" **
c  *************************************
      subroutine makeinit_cvb()
      implicit real*8 (a-h,o-z)
#include "make_cvb.fh"

      nobj=0
      ndep_ij=0
      ndep_ji=0
      ioffs(1)=0
      joffs(1)=0
      mustdeclare=.false.
      iprint=0
      return
      end
