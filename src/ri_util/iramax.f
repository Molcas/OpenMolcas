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
* Copyright (C) Francesco Aquilante                                    *
************************************************************************
      integer function irAmax(n,a,inc)
c
c     Finds the index of element having max. absolute value.
c
c     Author:  F. Aquilante
c
      implicit real*8 (a-h,o-z)
      REAL*8 a(*)
      integer n,inc
c
      iramax = 0
      if( n.lt.1 .or. inc.le.0 )return
      iramax = 1
      if(n.eq.1)return
      if(inc.eq.1)go to 20
c
c        code for increment not equal to 1
c
      ix = 1
      smax = abs(a(1))
      ix = ix + inc
      do 10 i = 2,n
         if(abs(a(ix)).le.smax) go to 5
         iramax = i
         smax = abs(a(ix))
    5    ix = ix + inc
   10 continue
      return
c
c        code for increment equal to 1
c
   20 smax = abs(a(1))
      do 30 i = 2,n
         if(abs(a(i)).le.smax) go to 30
         iramax = i
         smax = abs(a(i))
   30 continue
      return
      end
