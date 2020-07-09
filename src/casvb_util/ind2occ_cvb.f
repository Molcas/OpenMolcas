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
      integer function indget_cvb(iminor, idim, nel, ixmin)
      dimension ixmin(0:nel,0:idim),iminor(nel)
      indget_cvb=1
      iacc=0
      do 1100 i=1,nel
      if(iminor(i).eq.1)then
        iacc=iacc+1
        indget_cvb=indget_cvb+ixmin(i-1,iacc)
      endif
1100  continue
      return
      end
      integer function minind_cvb(iminor, idim, nel, ixmin)
      dimension ixmin(0:nel,0:idim),iminor(idim)
      minind_cvb=1
      do 1100 i=1,idim
      minind_cvb=minind_cvb+ixmin(iminor(i)-1,i)
1100  continue
      return
      end
