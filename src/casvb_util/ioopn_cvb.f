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
      subroutine ioopn_cvb(fn,lu)
      implicit real*8(a-h,o-z)
      character*(*) fn
#include "fio.fh"

c  Close first in case file is 'hanging' from previous session:
      if(isOpen(lu).ne.0)call daclos(lu)
      call daname_wa(lu,fn)
      end
