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
      subroutine weightfl_cvb(ix,n,nel)
      dimension ix(0:nel,0:n)
      call izero(ix,(n+1)*(nel+1))
      ix(0,0)=1
      do 1200 iel=1,nel
      do 1200 ik=max(iel-nel+n,0),min(iel,n)
      if(ik.ne.0)then
        ix(iel,ik)=ix(iel-1,ik) + ix(iel-1,ik-1)
      else
        ix(iel,ik)=ix(iel-1,ik)
      endif
1200  continue
      return
      end
