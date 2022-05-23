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
      subroutine loopstr0_cvb(iocc,index,nel,norb)
      implicit real*8 (a-h,o-z)
      dimension iocc(nel)

      if(nel.gt.norb)then
        write(6,*)' Error in loopstr0, nel > norb :',nel,norb
        call abend_cvb()
      endif
      index=1
      do 100 iel=1,nel
      iocc(iel)=iel
100   continue
      return
      end
