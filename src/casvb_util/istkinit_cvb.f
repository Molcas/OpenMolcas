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
c  ******************************
c  ** Simple stack of integers **
c  ******************************
      subroutine istkinit_cvb(iarr,n)
      implicit real*8(a-h,o-z)
      dimension iarr(n)

      if(n.lt.2)then
        write(6,*)' Too small dimension in ISTKINIT_CVB :',n
        call abend_cvb()
      endif
      iarr(1)=n
      iarr(2)=2
      return
      end
