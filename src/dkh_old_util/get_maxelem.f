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
* Copyright (C) 2004,2005, Alexander Wolf                              *
*               2004,2005, Markus Reiher                               *
************************************************************************
      subroutine get_maxelem (isize,array,pos,maxelem)
c
************************************************************************
c
c   This SR belongs to the test suite 'parser2'
c
c   written by:  Alexander Wolf and Markus Reiher  (Univ. Jena)
c
c   version:  2.0.1
c
c   last modified: 21.07.2005
c
c   first version: 10.04.2004  (Theoretical Chemistry, Univ. Bonn)
c
c************************************************************************
c
      implicit none
      integer isize,pos,i
      real*8 array(isize),maxelem
      intrinsic ABS
c
      pos=0
      maxelem=0.0d0
      do 10 i=1,isize
        if (abs(array(i)).gt.maxelem) then
          pos=i
          maxelem=abs(array(i))
        endif
  10  continue
c
      return
      end
