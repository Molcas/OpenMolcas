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
      subroutine svd_cvb(a,val,vec,vmat,n1,n2)
      implicit real*8 (a-h,o-z)
#include "malloc_cvb.fh"
      dimension a(n1,n2),val(n2),vec(n1,n2),vmat(n2,n2)

      n12=max(n1,n2)
      i1 = mstackr_cvb(n12*n2)
      i2 = mstackr_cvb(n2)
      i3 = mstackr_cvb(n12*n2)
      i4 = mstackr_cvb(n12*n2)
      i5 = mstackr_cvb(n2)
      i6 = mstacki_cvb(n2)
      call svd2_cvb(a,val,vec,vmat,n1,n2,n12,
     >  w(i1),w(i2),w(i3),w(i4),w(i5),iw(i6))
      call mfreer_cvb(i1)
      return
      end
