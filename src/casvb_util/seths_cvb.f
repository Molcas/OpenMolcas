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

      call seth_cvb(n,1)
      lenarr=len(arr(1))
      do 100 i=1,n
      do 100 j=1,lenarr
100   call seth_cvb(ichar(arr(i)(j:j)),1)
      return
      end

      subroutine geths_cvb(arr,n)
      implicit real*8(a-h,o-z)
      character*(*) arr(n)
      call geth_cvb(n,1)
      lenarr=len(arr(1))
      do 200 i=1,n
      do 200 j=1,lenarr
      call geth_cvb(iret,1)
200   arr(i)(j:j)=char(iret)
      return
      end
