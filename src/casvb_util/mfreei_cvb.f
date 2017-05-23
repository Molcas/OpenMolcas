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
      subroutine mfreei_cvb(ipoint)
      implicit real*8 (a-h,o-z)
#include "idbl_cvb.fh"
#include "memman_cvb.fh"

      if(memdebug)write(6,*)'   Enter mfreei: pointer :',ipoint
      iraddr=(ipoint-1)/idbl+1
      call mfreer_cvb(iraddr)
      return
      end
