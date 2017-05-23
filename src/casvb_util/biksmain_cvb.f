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
      subroutine biksmain_cvb(aikcof,bikcof,
     > nel,nalf,ndet,ifns,kbasis,share,iprint)
      implicit real*8 (a-h,o-z)
      logical share
      dimension aikcof(ndet,ifns),bikcof(ndet,ifns)
#include "malloc_cvb.fh"
      save one
      data one/1d0/

      if(nel.eq.0.and.kbasis.ne.6) then
        bikcof(1,1)=one
        aikcof(1,1)=one
        return
      endif

      nbet=nel-nalf
      np1=nel+1
      nalf1=nalf+1
      nbet2=nbet+nbet
      nswpdim=2**nbet

      i1 = mstacki_cvb(np1)
      i2 = mstacki_cvb(np1)
      i3 = mstacki_cvb(np1)
      i4 = mstacki_cvb(nbet2+1)
      i5 = mstacki_cvb(nbet2+1)
      i6 = mstacki_cvb(nbet2+1)
      i7 = mstacki_cvb(nbet*nswpdim)
      i8 = mstacki_cvb(nel)
      i9 = mstacki_cvb(nel)
      i10= mstacki_cvb(np1*nalf1)
      i11= mstacki_cvb(np1*nalf1)
      i12= mstacki_cvb(nel)
      i13= mstacki_cvb(nalf)
      i14= mstacki_cvb(nbet)
      call rumer_cvb(bikcof,
     >  nel,nalf,nbet,ndet,ifns,kbasis,iprint,nswpdim,
     >  iw(i1),iw(i2),iw(i3),
     >  iw(i4),iw(i5),iw(i6),iw(i7),iw(i8),iw(i9),
     >  iw(i10),iw(i11),iw(i12),
     >  iw(i13),iw(i14))
      call mfreei_cvb(i1)

      if(kbasis.eq.1.or.kbasis.eq.5)
     >  call kotani_cvb(bikcof,ndet,ifns)

      if(kbasis.eq.5)then
        i1 = mstacki_cvb(np1)
        i2 = mstacki_cvb(np1)
        i3 = mstacki_cvb(np1)
        i4 = mstacki_cvb(np1)
        i5 = mstacki_cvb(np1)
        i6 = mstacki_cvb(np1)
        i7 = mstacki_cvb(nbet2)
        i8 = mstacki_cvb(nbet2)
        i9 = mstacki_cvb(nel)
        i10= mstacki_cvb(nel)
        i11= mstacki_cvb(nalf)
        i12= mstacki_cvb(np1*nalf1)
        i13= mstacki_cvb(np1*nalf1)
        i14= mstacki_cvb(nel)
        i15= mstackr_cvb(ndet)
        call projspn_cvb(bikcof,
     >    nel,nalf,nbet,ndet,ifns,
     >    iw(i1),iw(i2),iw(i3),iw(i4),iw(i5),iw(i6),
     >    iw(i7),iw(i8),iw(i9),iw(i10),
     >    iw(i11),
     >    iw(i12),iw(i13),iw(i14),w(i15))
        call mfreei_cvb(i1)
      endif

      if(kbasis.eq.2)then
        i1 = mstacki_cvb(np1)
        i2 = mstacki_cvb(np1)
        i3 = mstacki_cvb(np1)
        i4 = mstacki_cvb(nel)
        i5 = mstacki_cvb(nel)
        i6 = mstacki_cvb(np1*nalf1)
        i7 = mstacki_cvb(nalf)
        i8 = mstacki_cvb(nbet)
        i9 = mstacki_cvb(ifns)
        call serber_cvb(bikcof,
     >    nel,nalf,nbet,ndet,ifns,
     >    iw(i1),iw(i2),iw(i3),iw(i4),iw(i5),iw(i6),
     >    iw(i7),iw(i8),iw(i9))
        call mfreei_cvb(i1)
      endif

      if(kbasis.gt.2.and.kbasis.ne.6)then
        i1 = mstackr_cvb(ifns*ifns)
      else
        i1 = mstackr_cvb(0)
      endif
      call aikcof_cvb(aikcof,bikcof,
     >  ndet,ifns,kbasis,share,
     >  w(i1))
      call mfreer_cvb(i1)

      return
      end
