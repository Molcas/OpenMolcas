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
      subroutine o7a_cvb(nparm)
      implicit real*8 (a-h,o-z)
#include "opt2_cvb.fh"
      save one
      data one/1d0/

      call ddnewopt_cvb()
      have_solved_it=.false.
      call ddguess_cvb([one],1,0)
      return
c Avoid unused argument warnings
      if (.false.) call Unused_integer(nparm)
      end
