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
c  ** Memory allocation and the like **
c  ************************************
      subroutine setidbl_cvb()
      implicit real*8(a-h,o-z)
#include "idbl_cvb.fh"
#include "SysDef.fh"

      idbl=RtoI
      return
      end
      integer function idbl_cvb(nreals)
      implicit real*8(a-h,o-z)
#include "idbl_cvb.fh"

      idbl_cvb=nreals*idbl
      return
      end
      integer function ihlf_cvb(nints)
      implicit real*8(a-h,o-z)
#include "idbl_cvb.fh"

      ihlf_cvb=(nints+idbl-1)/idbl
      return
      end
c
c  -- Initialization of casvb memory manager ---
c
