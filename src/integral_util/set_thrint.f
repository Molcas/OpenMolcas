************************************************************************
* This file is part of OpenMolcas.                                     *
*                                                                      *
* OpenMolcas is free software; you can redistribute it and/or modify   *
* it under the terms of the GNU Lesser General Public License, v. 2.1. *
* OpenMolcas is distributed in the hope that it will be useful, but it *
* is provided "as is" and without any express or implied warranties.   *
* For more details see the full text of the license in the file        *
* LICENSE or in <http://www.gnu.org/licenses/>.                        *
************************************************************************
************************************************************************
*                                                                      *
*  Subroutine set_thrint:     Puts integral thresholds into info.fh    *
*                                                                      *
************************************************************************
      subroutine set_thrint(thr,cut)
c-----------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
#include "itmax.fh"
#include "info.fh"
      thrint=max(1.d-14,thr)
      cutint=max(1.d-14,cut)
      return
      end
