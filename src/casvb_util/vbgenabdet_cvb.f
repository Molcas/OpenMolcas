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
      subroutine vbgenabdet_cvb(idetavb,idetbvb,
     >  iconfs,nconf,nconfion,
     >  ndetvb,nel,noe,
     >  nalf,nbet,norb)
      implicit real*8 (a-h,o-z)
#include "malloc_cvb.fh"
      dimension idetavb(ndetvb),idetbvb(ndetvb)
      dimension iconfs(noe,nconf),nconfion(0:nel)

      i1 = mstacki_cvb((norb+1)*(nalf+1))
      i2 = mstacki_cvb((norb+1)*(nbet+1))
      i3 = mstacki_cvb(norb+1)
      i4 = mstacki_cvb(norb+1)
      i5 = mstacki_cvb(norb)
      i6 = mstacki_cvb(norb)
      call vbgenabdet2_cvb(idetavb,idetbvb,
     >  iconfs,nconf,nconfion,
     >  ndetvb,nel,noe,
     >  nalf,nbet,norb,
     >  iw(i1),iw(i2),iw(i3),iw(i4),
     >  iw(i5),iw(i6))
      call mfreei_cvb(i1)
      return
      end
