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
      subroutine loopstr_cvb(iocc,index,nel,norb)
      implicit real*8 (a-h,o-z)
      dimension iocc(nel)

      index=index+1
c  Find electron for which orbital number can be increased :
      do 100 iel=1,nel-1
      if(iocc(iel+1).gt.iocc(iel)+1)goto 200
100   continue
      iel=nel
      if(iocc(iel).lt.norb)goto 200
      call loopstr0_cvb(iocc,index,nel,norb)
      return
200   iocc(iel)=iocc(iel)+1
      do 300 jel=1,iel-1
      iocc(jel)=jel
300   continue
      return
      end
