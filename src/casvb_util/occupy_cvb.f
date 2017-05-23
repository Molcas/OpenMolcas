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
      subroutine occupy_cvb(nk,nel,locc,lunocc)
      dimension nk(0:nel),locc(*),lunocc(*)

      iocc=0
      iunocc=0
      do 2100 iel=1,nel
      if(nk(iel)-nk(iel-1) .eq. 1)then
        iocc=iocc+1
        locc(iocc)=iel
      elseif(nk(iel)-nk(iel-1) .eq. 0)then
        iunocc=iunocc+1
        lunocc(iunocc)=iel
      else
        write(6,*)' Error in graphical indexing routine!'
        call abend_cvb()
      endif
2100  continue
      return
      end
