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
* Copyright (C) 1996-2006, Thorstein Thorsteinsson                     *
*               1996-2006, David L. Cooper                             *
************************************************************************
      subroutine sminus_cvb(bikfrom,bikto,
     >  nel,nalffrom,nalfto,nvec)
      implicit real*8 (a-h,o-z)
      dimension bikfrom(*),bikto(*)
#include "malloc_cvb.fh"

      call ab2asc_cvb(bikfrom,nvec,nel,nalffrom)

      do 100 ialfto=nalffrom-1,nalfto,-1
      ialffrom=ialfto+1
      i1 = mstacki_cvb((nel+1)*(ialffrom))
      i2 = mstacki_cvb(ialffrom)
      i3 = mstacki_cvb(ialfto)
      ndetfrom=ndet_cvb(nel,ialffrom)
      ndetto=ndet_cvb(nel,ialfto)
      if(nalffrom.eq.nalfto+1)then
        call sminus2_cvb(bikfrom,bikto,
     >    nel,ialffrom,ndetfrom,ialfto,ndetto,nvec,
     >    iw(i1),iw(i2),iw(i3))
      elseif(ialfto.eq.nalffrom-1)then
        i4 = mheapr_cvb(ndetto*nvec)
        call sminus2_cvb(bikfrom,w(i4),
     >    nel,ialffrom,ndetfrom,ialfto,ndetto,nvec,
     >    iw(i1),iw(i2),iw(i3))
      elseif(ialfto.eq.nalfto)then
        call sminus2_cvb(w(i4),bikto,
     >    nel,ialffrom,ndetfrom,ialfto,ndetto,nvec,
     >    iw(i1),iw(i2),iw(i3))
        call mhpfreer_cvb(i4)
      else
        i5 = mheapr_cvb(ndetto*nvec)
        call sminus2_cvb(w(i4),w(i5),
     >    nel,ialffrom,ndetfrom,ialfto,ndetto,nvec,
     >    iw(i1),iw(i2),iw(i3))
        call mhpfreer_cvb(i4)
        i4=i5
      endif
      call mfreei_cvb(i1)
100   continue

      call asc2ab_cvb(bikto,nvec,nel,nalfto)
c  Now try to retain normalization ...
      ndetfrom=ndet_cvb(nel,nalffrom)
      ndetto=ndet_cvb(nel,nalfto)
      do 200 ivec=1,nvec
      cnrmfrom=dnrm2_(ndetfrom,bikfrom(1+(ivec-1)*ndetfrom),1)
      cnrmto=dnrm2_(ndetto,bikto(1+(ivec-1)*ndetto),1)
      if(cnrmto.gt.1d-10)
     >  call dscal_(ndetto,cnrmfrom/cnrmto,bikto(1+(ivec-1)*ndetto),1)
200   continue
      return
      end
