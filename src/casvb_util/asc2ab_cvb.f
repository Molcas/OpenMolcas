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
      subroutine asc2ab_cvb(detvec,nvec,nel,nalf)
      implicit real*8 (a-h,o-z)
#include "malloc_cvb.fh"
      dimension detvec(*)


      call icomb_cvb(nel,nalf,ndet)
      nbet=nel-nalf
      i1 = mstacki_cvb(nel+1)
      i2 = mstacki_cvb(nel+1)
      i3 = mstacki_cvb(nel+1)
      i4 = mstacki_cvb((nel+1)*(nalf+1))
      i5 = mstacki_cvb(nel)
      call asc2ab2_cvb(detvec,nvec,nel,nalf,
     >  nbet,ndet,
     >  iw(i1),iw(i2),iw(i3),iw(i4),iw(i5))
      call mfreei_cvb(i1)
      return
      end
      subroutine ab2asc_cvb(detvec,nvec,nel,nalf)
      implicit real*8 (a-h,o-z)
#include "malloc_cvb.fh"
      dimension detvec(*)


      call icomb_cvb(nel,nalf,ndet)
      nbet=nel-nalf
      i1 = mstacki_cvb(nel+1)
      i2 = mstacki_cvb(nel+1)
      i3 = mstacki_cvb(nel+1)
      i4 = mstacki_cvb((nel+1)*(nalf+1))
      i5 = mstacki_cvb(nel)
      call asc2ab2_cvb(detvec,nvec,nel,nalf,
     >  nbet,ndet,
     >  iw(i1),iw(i2),iw(i3),iw(i4),iw(i5))
      call mfreei_cvb(i1)
      return
      end
