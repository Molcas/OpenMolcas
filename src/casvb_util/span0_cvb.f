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
      subroutine span0_cvb(nvecmx1,n)
      implicit real*8 (a-h,o-z)
      common /span_comcvb/iaddr,nvecmx,nvtot
      save nmult
      data nmult/5/

      nvecmx=min(nmult*nvecmx1,mavailr_cvb()/n)
      if(nvecmx.le.0)then
        write(6,*)' Not enough vectors in SPAN0_CVB!',nvecmx
        write(6,*)' Remaining memory :',mavailr_cvb()
        write(6,*)' Max number of vectors :',nvecmx1
        call abend_cvb()
      endif
      iaddr = mstackr_cvb(n*nvecmx)
      nvtot=0
      return
      end
