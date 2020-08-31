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
c  **************************************
c  ** General-purpose utility routines **
c  **************************************
c  ********************************************
c  ** Simple vector/matrix-handling routines **
c  ********************************************
      subroutine mxfold_cvb(a,avec,n)
      implicit real*8 (a-h,o-z)
      dimension a(n,n),avec(n*(n-1))
      iprm=0
      do 100 i=1,n
      do 101 j=1,n
      if(j.ne.i)then
        iprm=iprm+1
        avec(iprm)=a(j,i)
      endif
101   continue
100   continue
      return
      end
