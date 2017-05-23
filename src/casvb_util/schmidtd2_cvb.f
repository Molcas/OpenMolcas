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
      subroutine schmidtd2_cvb(c1,sxc1,nvec1,c2,nvec2,n)
      implicit real*8 (a-h,o-z)
      dimension c1(n,nvec1),sxc1(n,nvec1),c2(n,nvec2)

      do 100 i=1,nvec2
      do 100 j=1,nvec1
100   call daxpy_(n,-ddot_(n,c2(1,i),1,sxc1(1,j),1),c1(1,j),1,c2(1,i),1)
      return
      end
