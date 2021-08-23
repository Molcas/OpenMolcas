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
c  *********************************************************************
c  *                                                                   *
c  *  CIWEIGHT  := Chirgwin-Couson/Lowdin/inverse-overlap weights of   *
c  *            := full CASSCF vector and residual.                    *
c  *                                                                   *
c  *********************************************************************
      subroutine ciweight_cvb(civec,civbs,civb,citmp,vec5,
     >  orbs,sorbs,orbinv,owrk,gjorb,gjorb2,gjorb3)
      implicit real*8 (a-h,o-z)
#include "main_cvb.fh"
#include "optze_cvb.fh"
#include "files_cvb.fh"
#include "print_cvb.fh"

#include "WrkSpc.fh"
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
      ncnfcas=ncnfcas+iretval1*iretval2
100   continue

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
     >  work(iaddr_ci(icitmp)),work(iaddr_ci(icivbs)),
     >  work(iaddr_ci(icivec)),
     >  work(iaddr_ci(icivb)),work(iaddr_ci(ivec5)),
     >  work(k1),work(k2),work(k3),work(k4),work(k5),work(k6),
     >  iwork(k7),iwork(k8),iwork(k9),iwork(k10),iwork(k11),iwork(k12),
     >  iwork(k13),iwork(k14),iwork(k15),iwork(k16),iwork(k17),
     >  iwork(k18),
     >  iwork(k19),iwork(k20),iwork(k21),iwork(k22),iwork(k23),
     >  iwork(k24),
     >  iwork(k25),iwork(k26),iwork(k27),iwork(k28),iwork(k29),
     >  iwork(k30),
     >  work(k31),work(k32),iwork(k33),iwork(k34),
     >  ionmin,ionmax,mxrem,mxsng,mxasg,ncnfcas,mxdetcas)
      call mfreer_cvb(k1)
      return
      end
