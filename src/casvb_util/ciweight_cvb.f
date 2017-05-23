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
c  *********************************************************************
c  *                                                                   *
c  *  CIWEIGHT  := Chirgwin-Couson/Lowdin/inverse-overlap weights of   *
c  *            := full CASSCF vector and residual.                    *
c  *                                                                   *
c  *********************************************************************
      subroutine ciweight_cvb(civec,civbs,civb,citmp,vec5,
     >  orbs,sorbs,orbinv,owrk,gjorb,gjorb2,gjorb3)
      implicit real*8 (a-h,o-z)
#include "ext_cvb.fh"
#include "main_cvb.fh"
#include "optze_cvb.fh"
#include "files_cvb.fh"
#include "print_cvb.fh"

#include "malloc_cvb.fh"
      dimension orbs(norb,norb),sorbs(norb,norb)
      dimension orbinv(norb,norb),owrk(norb,norb)
      dimension gjorb(*),gjorb2(*),gjorb3(*),civec(*),civbs(*),
     > civb(*),citmp(*),vec5(*)

      ionmin=max(nel-norb,0)
      ionmax=nbet
      mxrem=norb-ionmin
      mxsng=nel-2*ionmin
      mxasg=nalf-ionmin
      call icomb_cvb(mxsng,mxasg,mxdetcas)
c  Work out number of configurations in CASSCF vector :
      ncnfcas=0
      do 100 ion=ionmin,ionmax
      call icomb_cvb(norb,ion,iretval1)
      call icomb_cvb(norb-ion,nel-2*ion,iretval2)
100   ncnfcas=ncnfcas+iretval1*iretval2

      k1 = mstackr_cvb(ionmax-ionmin+1)
      k2 = mstackr_cvb(ionmax-ionmin+1)
      k3 = mstackr_cvb(ionmax-ionmin+1)
      k4 = mstackr_cvb(ionmax-ionmin+1)
      k5 = mstackr_cvb(ionmax-ionmin+1)
      k6 = mstackr_cvb(ionmax-ionmin+1)
      k7 = mstacki_cvb(norb+1)
      k8 = mstacki_cvb(norb+1)
      k9 = mstacki_cvb((norb+1)*(nalf+1))
      k10= mstacki_cvb((norb+1)*(nbet+1))
      k11= mstacki_cvb(norb)
      k12= mstacki_cvb(norb)
      k13= mstacki_cvb(norb+1)
      k14= mstacki_cvb(norb+1)
      k15= mstacki_cvb(norb+1)
      k16= mstacki_cvb((norb+1)*(ionmax+1))
      k17= mstacki_cvb(norb)
      k18= mstacki_cvb(norb)
      k19= mstacki_cvb(norb+1)
      k20= mstacki_cvb(norb+1)
      k21= mstacki_cvb(norb+1)
      k22= mstacki_cvb((mxrem+1)*(mxsng+1))
      k23= mstacki_cvb(norb)
      k24= mstacki_cvb(norb)
      k25= mstacki_cvb(norb+1)
      k26= mstacki_cvb(norb+1)
      k27= mstacki_cvb(norb+1)
      k28= mstacki_cvb((mxsng+1)*(mxasg+1))
      k29= mstacki_cvb(norb)
      k30= mstacki_cvb(norb)
      if(mod(iciweights,8).gt.3)then
        k31= mstackr_cvb(ncnfcas)
        k32= mstackr_cvb(ncnfcas)
        k33= mstacki_cvb(mxdetcas)
        k34= mstacki_cvb(mxdetcas)
      else
        k31= 0
        k32= 0
        k33= 0
        k34= 0
      endif
      icivec=nint(civec(1))
      icivbs=nint(civbs(1))
      icivb=nint(civb(1))
      icitmp=nint(citmp(1))
      ivec5=nint(vec5(1))
      call ciweight2_cvb(civec,civbs,civb,citmp,vec5,
     >  orbs,sorbs,orbinv,owrk,gjorb,gjorb2,gjorb3,
     >  w(iaddr_ci(icitmp)),w(iaddr_ci(icivbs)),w(iaddr_ci(icivec)),
     >  w(iaddr_ci(icivb)),w(iaddr_ci(ivec5)),
     >  w(k1),w(k2),w(k3),w(k4),w(k5),w(k6),
     >  iw(k7),iw(k8),iw(k9),iw(k10),iw(k11),iw(k12),
     >  iw(k13),iw(k14),iw(k15),iw(k16),iw(k17),iw(k18),
     >  iw(k19),iw(k20),iw(k21),iw(k22),iw(k23),iw(k24),
     >  iw(k25),iw(k26),iw(k27),iw(k28),iw(k29),iw(k30),
     >  w(k31),w(k32),iw(k33),iw(k34),
     >  ionmin,ionmax,mxrem,mxsng,mxasg,ncnfcas,mxdetcas)
      call mfreer_cvb(k1)
      return
      end
