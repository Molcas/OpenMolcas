************************************************************************
* This file is part of OpenMolcas.                                     *
*                                                                      *
* OpenMolcas is free software; you can redistribute it and/or modify   *
* it under the terms of the GNU Lesser General Public License, v. 2.1. *
* OpenMolcas is distributed in the hope that it will be useful, but it *
* is provided "as is" and without any express or implied warranties.   *
* For more details see the full text of the license in the file        *
* LICENSE or in <http://www.gnu.org/licenses/>.                        *
************************************************************************
      subroutine import_cho(numcho_pt2,infvec_n2_pt2,maxvec_pt2)
      implicit none
#include "cholesky.fh"
      integer numcho_pt2(8), infvec_n2_pt2, maxvec_pt2
      integer i
      do i=1,nsym
        numcho_pt2(i)=numcho(i)
      end do
      do i=nsym+1,8
        numcho_pt2(i)=0
      end do
      maxvec_pt2=maxvec
      infvec_N2_pt2=infvec_N2
      return
      end
