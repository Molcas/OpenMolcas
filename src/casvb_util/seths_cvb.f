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
      subroutine seths_cvb(arr,n)
      implicit real*8(a-h,o-z)
      character*(*) arr(n)

      call seth_cvb([n],1)
      lenarr=len(arr(1))
      do 100 i=1,n
      do 101 j=1,lenarr
      call seth_cvb([ichar(arr(i)(j:j))],1)
101   continue
100   continue
      return
      end

      subroutine geths_cvb(arr,n)
      implicit real*8(a-h,o-z)
      character*(*) arr(n)
      dimension iaux(1)
      call geth_cvb(iaux,1)
      n=iaux(1)
      lenarr=len(arr(1))
      do 200 i=1,n
      do 201 j=1,lenarr
      call geth_cvb(iaux,1)
      iret=iaux(1)
      arr(i)(j:j)=char(iret)
201   continue
200   continue
      return
      end
